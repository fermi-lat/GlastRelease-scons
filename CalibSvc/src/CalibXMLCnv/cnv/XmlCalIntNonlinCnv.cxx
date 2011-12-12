// $Header$

#include <string>
#include "XmlCalIntNonlinCnv.h"

#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GenericAddress.h"

#include "CalibSvc/ICalibXmlSvc.h"
#include "CalibSvc/ICalibMetaCnvSvc.h"

#include "CalibData/Cal/CalCalibIntNonlin.h"

// Probably need the following:
#include "CalibData/DacCol.h"

#include "CalibData/CalibTime.h"
#include "xmlBase/Dom.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

//static CnvFactory<XmlCalIntNonlinCnv> s_factory;
//const  ICnvFactory& XmlCalIntNonlinCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(XmlCalIntNonlinCnv);



XmlCalIntNonlinCnv::XmlCalIntNonlinCnv( ISvcLocator* svc) :
  XmlCalBaseCnv(svc, CLID_Calib_CAL_IntNonlin) { 
}


const CLID& XmlCalIntNonlinCnv::objType() const {
  return CLID_Calib_CAL_IntNonlin;
}

const CLID& XmlCalIntNonlinCnv::classID() {
  return CLID_Calib_CAL_IntNonlin;
}

namespace {
  /// Local utility which knows how to get the information out of a
  /// <calLightAsym> element and make a CalibData::IntNonlin with it
  CalibData::IntNonlin* processRange(DOMElement* intNonlinElt) {
    using xmlBase::Dom;

    // Could check here to make sure it really is a <intNonlin>

    float err;
    std::vector<float> vals;
    std::vector<float> sdacs;
    unsigned nSdacs = 0;
    try {
      xmlBase::Dom::getFloatsAttribute(intNonlinElt, "values", vals);
      nSdacs =  xmlBase::Dom::getFloatsAttribute(intNonlinElt, "sdacs", sdacs);
      err = xmlBase::Dom::getDoubleAttribute(intNonlinElt, "error");
    }
    catch (xmlBase::DomException& ex) {
      std::cerr << "From CalibSvc::XmlCalIntNonlinCnv::processRange" 
                << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      std::cerr.flush();
      return 0;
    }
    if (nSdacs > 0) {
      if (nSdacs == vals.size()) {
        return new CalibData::IntNonlin(&vals, err, &sdacs);
      }
      else {
        std::cerr << "From CalibSvc::XmlCalIntNonlinCnv::processRange" 
                  << std::endl;
        std::cerr << "#dacs (" << nSdacs << ") != #values (" 
                  << vals.size() << ")" << std::endl;
        std::cerr.flush();
        return 0;
      }
    }
    else return new CalibData::IntNonlin(&vals, err, 0);
  }
}

// Create our specific object
StatusCode XmlCalIntNonlinCnv::i_createObj(const DOMElement* docElt, 
                                     DataObject*& refpObject)
{
  using xmlBase::Dom;
  using CalibData::IntNonlin;
  using CalibData::DacCol;

  unsigned nRow, nCol, nLayer, nXtal, nFace, nRange, nDacCol;
  // need dimensions to call the constructor
  StatusCode status = 
    readDimension(docElt, nRow, nCol, nLayer, nXtal, nFace, nRange,
                  &nDacCol);
  if (status == StatusCode::FAILURE) return status;
  
  // refpObject
  CalibData::CalCalibIntNonlin* pObj = 
    new CalibData::CalCalibIntNonlin(nRow, nCol, nLayer, nXtal, nFace, nRange,
                                     nDacCol);
  refpObject = pObj;
  if (!pObj) return StatusCode::FAILURE;

  setBaseInfo(pObj);

  DOMElement* rangeElt = findFirstRange(docElt);

  while (rangeElt != 0 ) {
    IntNonlin* pIntNonlin = processRange(rangeElt);
    if (pIntNonlin) {
      pObj->putRange(m_nRow, m_nCol, m_nLayer, m_nXtal, m_nRange, 
                     m_nFace, pIntNonlin);
      delete pIntNonlin;
      rangeElt = findNextRange(rangeElt);
    }
    else {
      std::cerr << 
        "Bad specification for (towerRow,towerCol,layer,xtal,range,face)"
                << std::endl << "(" << m_nRow 
                << "," << m_nCol << "," << m_nLayer << "," << m_nXtal 
                << "," << m_nRange << "," << m_nFace << ")" << std::endl;
      std::cerr << "Exiting.." << std::endl;
      std::cerr.flush();
      exit(1);
    }
  }

  // Also have to handle dac collections        <<<<-------
  DOMElement* dacColElt = findFirstDacCol(docElt);

  while (dacColElt != 0 ) {
    unsigned range;
    DacCol* pDacCol = processDacCol(dacColElt, &range);
    pObj->putDacCol(range, pDacCol);
    delete pDacCol;

    dacColElt = findNextDacCol(dacColElt);
  }
  return StatusCode::SUCCESS;
}



#ifndef CALIBSVC_XmlAcdCnv_h
#define CALIBSVC_XmlAcdCnv_h

#include "XmlAcdBaseCnv.h"

#include "xmlBase/Dom.h"
#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/Acd/AcdCalib.h"

namespace CalibData {
  class AcdCalibDescription;
}

template <class T>
class XmlAcdCnv : public XmlAcdBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory< XmlAcdCnv<T> >;

public:
  const CLID& objType() const {
    return T::calibCLID();
  }
  static const CLID& classID() {
    return T::calibCLID();
  }

protected:

  XmlAcdCnv(ISvcLocator* svc)
    :XmlAcdBaseCnv(svc,T::calibCLID()){
  }

public:
  virtual ~XmlAcdCnv() {}       // most likely nothing to do 

  T* processPmt(DOMElement* pmtElt, const CalibData::AcdCalibDescription& desc) {
    /// Local utility which knows how to get the information out of a
    /// <acdXXX> element and make a CalibData::AcdXXX with it
    using xmlBase::Dom;
    
    // Element we're interested in is child of <pmt>
    DOMElement* calibElt = xmlBase::Dom::getFirstChildElement(pmtElt);

    // Could check here to make sure it really is an <acdPed>
    int status;
    std::vector<float> vals;
    try {    
      status = xmlBase::Dom::getIntAttribute(calibElt, "status");
      for ( int i(0); i < desc.nVar(); i++ ) {
	float varVal = xmlBase::Dom::getDoubleAttribute(calibElt, desc.getVarName(i).c_str());
	vals.push_back(varVal);
      }
    }
    catch (xmlBase::DomException ex) {
      std::cerr << "From CalibSvc::XmlAcdCnv::processPmt" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }
    return new T(desc,vals,(CalibData::AcdCalibObj::STATUS)status);
  }
  
  virtual StatusCode i_createObj(const DOMElement* element,
                                 DataObject*& refpObject) {
    
    using xmlBase::Dom;
    
    unsigned nFace, nRow, nCol, nPmt, nDet, nNA;
    // need dimensions to call the constructor
    StatusCode status = 
      readAcdDimension(element, nDet, nFace, nRow, nCol, nNA, nPmt);
    if (status == StatusCode::FAILURE) return status;
    
    // refpObject
    const CalibData::AcdCalibDescription* desc = CalibData::AcdCalibDescription::getDesc( T::calibType(), -1 );
    
    CalibData::AcdCalib<T>* pObj = 
      new CalibData::AcdCalib<T>(*desc, nFace, nRow, nCol, nNA, nPmt);
    refpObject = pObj;
    if (!pObj) return StatusCode::FAILURE;
    
    setBaseInfo(pObj);
    
    DOMElement* pmtElt = findFirstPmt(element);
    
    while (pmtElt != 0 ) {
      T* pPmt = processPmt(pmtElt,*desc);
      pObj->putPmt(m_id, m_nPmt, *pPmt);
      pmtElt = findNextPmt(pmtElt);
    }
    return StatusCode::SUCCESS;
  } 

};

#endif

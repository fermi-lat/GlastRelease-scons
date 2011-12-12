//  $Header$

#include <string>
#include "RootCalGainCnv.h"
#include "TTree.h"
#include "TFile.h"

#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GenericAddress.h"

#include "CalibSvc/ICalibRootSvc.h"    //maybe
#include "CalibSvc/ICalibMetaCnvSvc.h"

#include "CalibData/Cal/CalCalibGain.h"
// #include "CalibData/CalibTime.h"
#include "commonRootData/idents/CalXtalId.h"
#include "calibRootData/Cal/CalGainCol.h"
#include "idents/CalXtalId.h"

// Temporary.  Hope to find a better way to do this
#include "CalibData/CalibModel.h"

//static CnvFactory<RootCalGainCnv> s_factory;
//const  ICnvFactory& RootCalGainCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY(RootCalGainCnv);



RootCalGainCnv::RootCalGainCnv( ISvcLocator* svc) :
  RootCalBaseCnv(svc, CLID_Calib_CAL_ElecGain) { 
}


const CLID& RootCalGainCnv::objType() const {
  return CLID_Calib_CAL_ElecGain;
}

const CLID& RootCalGainCnv::classID() {
  return CLID_Calib_CAL_ElecGain;
}


StatusCode RootCalGainCnv::createRoot(const std::string& fname, 
                                      CalibData::CalibBase* pTDSObj) {
  // Make ROOT object of proper type
  //  calibRootData::CalGainCol gainCol;
  calibRootData::CalGainCol*  pGainCol = new calibRootData::CalGainCol;

  // This apparently is necessary.  Compiler doesn't like call to
  // RootBaseCnv::openWrite without it
  TObject* pTObj = pGainCol;

  std::string className("calibRootData::CalGainCol");

  // Open the file, create the branch
  StatusCode sc = openWrite(fname, className, pTObj);

  if (sc.isSuccess()) sc = fillRoot(pTDSObj, pGainCol);

  closeWrite();

  delete pGainCol;
  return sc;

}

StatusCode RootCalGainCnv::fillRoot(CalibData::CalibBase* pTDSObj, 
                                    TObject* pRootObj) {

  StatusCode sc;
  MsgStream log(msgSvc(), "RootCalGainCnv");

  // Make sure we've go sensible arguments
  calibRootData::CalGainCol* pRootGain = 
    dynamic_cast<calibRootData::CalGainCol* >(pRootObj);

  CalibData::CalCalibGain* pTDSGain = 
    dynamic_cast<CalibData::CalCalibGain* > (pTDSObj);

  if (!pRootGain) {
    log << "Cannot convert: Root object of wrong type " << endreq;
    return StatusCode::FAILURE;
  }

  if (!pTDSGain) {
    log << "Cannot convert: Calib TDS object of wrong type " << endreq;
    return StatusCode::FAILURE;
  }

  // Let base classes do their part
  sc = RootCalBaseCnv::fillRoot(pTDSObj, pRootGain->getDimension() );

  if (sc.isSuccess()) {
    // do our part. First store dimension information conveniently
    readDimension(pTDSGain);

    // For each allowed CalXtalId, see if TDS has a gain
    // object for it.  If so, make the Root version and append.
    // To avoid writing this horrible nested thing more than once,
    // should probably keep context of where we are in the loop
    // in RootCalBaseCnv and define a 'getNext' method there. 
    for (unsigned iRow = 0; iRow < m_nRow; iRow++) {
      for (unsigned iCol = 0; iCol < m_nCol; iCol++) {
        for (unsigned iLayer = 0; iLayer < m_nLayer; iLayer++) {
          for (unsigned iXtal = 0; iXtal < m_nXtal; iXtal++) {
            for (unsigned iFace = 0; iFace < m_nFace; iFace++) {
              for (unsigned iRange = 0; iRange < m_nRange; iRange++) {
                
                CalibData::RangeBase* pData = 
                  pTDSGain->getRange(iRow, iCol, iLayer, iXtal, iRange, iFace);
                if (pData) {
                  CalibData::Gain* pGain = 
                    dynamic_cast<CalibData::Gain*>(pData);
                  // had better not fail!
                  if (!pGain) continue;
                  CalXtalId id(m_nCol*iRow+iCol, iLayer, iXtal, iFace, iRange);
                  pRootGain->addChannel(id, pGain->getGain(),
                                        pGain->getSig());
                }
              }
            }
          }
        }
      }
    }

  }

  return sc;
}

StatusCode RootCalGainCnv::i_createObj(const std::string& fname,
                                       DataObject*& refpObject) {

  MsgStream log(msgSvc(), "RootCalGainCnv");
  calibRootData::CalGainCol* pCol = new calibRootData::CalGainCol;
  TObject* pTObj = pCol;
  //  StatusCode sc = openRead(fname, "calibRootData::CalGainCol", pTObj);
  // Note: RootBaseCnv::openRead now just opens the file, doesn't read in 
  // an "event" from the "Calib" tree.  
  StatusCode sc = openRead(fname);
  if (!sc.isSuccess() ) {
    delete pCol;
    return sc;
  }
  
  // Read in our object
  sc = readRootObj("calibRootData::CalGainCol", pTObj);
  // Must have dimensions even before we call constructor for TDS obj
  calibRootData::CalDimension* pDim = pCol->getDimension();

  if (!pDim) {
    delete pCol;
    return StatusCode::FAILURE;
  }
  CalibData::CalCalibGain* pObj = 
    new CalibData::CalCalibGain(pDim->getNRow(), pDim->getNCol(), 
                                pDim->getNLayer(), pDim->getNXtal(),
                                pDim->getNFace(), pDim->getNRange());
  if (!pObj) {
    delete pCol;
    return StatusCode::FAILURE;
  }
  refpObject = pObj;

  setBaseInfo(pObj);

  std::vector<calibRootData::CalGain>& rootGains = pCol->getGains();
  unsigned nGains = rootGains.size();
  log << MSG::DEBUG << "Copying gain information for " << nGains
      << " channels from ROOT to calibration TDS " << endreq;
  for (unsigned i = 0; i < nGains; i++) {
    CalibData::Gain* pGain = 
      new CalibData::Gain(rootGains[i].getGain(), rootGains[i].getSig());
    CalXtalId id = rootGains[i].getId();
    idents::CalXtalId identsId(id.getTower(), id.getLayer(), id.getColumn());
    pObj->putRange(identsId, id.getRange(), id.getFace(), pGain);
  }
  // All done with ROOT object.
  delete pCol;
  return StatusCode::SUCCESS;
}


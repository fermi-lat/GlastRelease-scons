// $Header$
#include "GaudiKernel/MsgStream.h"
#include "RootCalBaseCnv.h"
#include "calibRootData/Cal/CalDimension.h"
#include "CalibData/Cal/CalCalibBase.h"
#include "CalibSvc/IInstrumentName.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TObject.h"


RootCalBaseCnv::RootCalBaseCnv(ISvcLocator* svc, const CLID& clid) :
  RootBaseCnv(svc, clid), m_nRow(10000), m_nCol(10000), 
  m_nLayer(10000), m_nXtal(10000),
  m_nFace(10000), m_nRange(10000) {}

/// Store dimension information from TDS into our data members
StatusCode RootCalBaseCnv::readDimension(CalibData::CalCalibBase*  pCalBase) {
  m_nRow = pCalBase->getNTowerRow();
  m_nCol = pCalBase->getNTowerCol();
  m_nLayer = pCalBase->getNLayer();
  m_nXtal = pCalBase->getNXtal();
  m_nFace = pCalBase->getNFace();
  m_nRange = pCalBase->getNRange();
  m_nDacCol = pCalBase->getNDacCol();

  return StatusCode::SUCCESS;
}

StatusCode RootCalBaseCnv::fillRoot(CalibData::CalibBase* pTDSObj, 
                                    TObject* pRootObj) {

  MsgStream log(msgSvc(), "RootCalBaseCnv");

  // Let base class do its part
  //  StatusCode sc =   RootBaseCnv::fillRoot(pTDSObj, pRootObj);
  //  if (!(sc.isSuccess()) ) return sc;

  calibRootData::CalDimension* pCalDimension = 
    dynamic_cast<calibRootData::CalDimension* >(pRootObj);

  CalibData::CalCalibBase* pTDSCal = 
    dynamic_cast<CalibData::CalCalibBase* > (pTDSObj);

  if (!pCalDimension) {
    log << "Cannot convert: Root object of wrong type " << endreq;
    return StatusCode::FAILURE;
  }

  if (!pTDSCal) {
    log << "Cannot convert: Calib TDS object of wrong type " << endreq;
    return StatusCode::FAILURE;
  }
      
  pCalDimension->initialize(pTDSCal->getNTowerRow(), pTDSCal->getNTowerCol(),
                       pTDSCal->getNLayer(), pTDSCal->getNXtal(), 
                       pTDSCal->getNFace(), pTDSCal->getNRange(), 
                       pTDSCal->getNDacCol() );


  // Get instrument name from InstrumentName service
  TString instr = TString((m_instrSvc->getInstrumentName()).c_str());

  pCalDimension->setInstrument(instr);
  return StatusCode::SUCCESS;
}

/* TODO
   Change implementation to invoke RootBaseCnv::readRootObj
*/
 StatusCode RootCalBaseCnv::readRootObj(const std::string& branch, 
                                        TObject*& pCalib) {
  TTree* pTree = (TTree*)m_inFile->Get("Calib");

  pTree->SetBranchAddress(branch.c_str(), &pCalib);
  pTree->GetEvent();
  return StatusCode::SUCCESS;
 }


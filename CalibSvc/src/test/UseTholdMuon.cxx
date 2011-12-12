//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibSvc/ICalibPathSvc.h"
#include "CalibData/Cal/CalTholdMuon.h"
#include "CalibData/Cal/Xpos.h"
#include "idents/CalXtalId.h"   

/**
   @file UseTholdMuon.cxx
   Simple algorithm to test functioning of "the other" TDS, Cal new 
   (Sept. 04) tholdMuon calibration data.
*/

  /** 
   @class UseTholdMuon

   Algorithm exemplifying retrieval and use of CAL_TholdMuon calibration.
*/
class UseTholdMuon : public Algorithm {

public:
  UseTholdMuon(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper function called by execute
  void processNew(CalibData::CalTholdMuonCol* pNew, const std::string& path);

  void printValSigs(MsgStream* log, const std::vector<CalibData::ValSig>* v, 
                    std::string title);

  void printValSig(MsgStream* log, const CalibData::ValSig* v, 
                   std::string title);
  IDataProviderSvc* m_pCalibDataSvc;
  ICalibPathSvc*    m_pCalibPathSvc;
  std::string       m_path;
  int               m_ser;
};


/// Instantiation of a static factory to create instances of this algorithm
//static const AlgFactory<UseTholdMuon> Factory;
//const IAlgFactory& UseTholdMuonFactory = Factory;
DECLARE_ALGORITHM_FACTORY(UseTholdMuon);


UseTholdMuon::UseTholdMuon(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
  m_ser(-1)
{
  // Declare properties here.

}


StatusCode UseTholdMuon::initialize() {
  StatusCode sc;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initialize()" << endreq;

  // So far don't have any properties, but in case we do some day..
  setProperties();


  sc = service("CalibDataSvc", m_pCalibDataSvc, true);

  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Could not get IDataProviderSvc interface of CalibDataSvc" 
	<< endreq;
    return sc;
  }

  sc = service("CalibDataSvc", m_pCalibPathSvc, true);

  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Could not get ICalibPathSvc interface of CalibDataSvc" 
	<< endreq;
    return sc;
  }

  m_path = 
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_CAL_TholdMuon,
                                  std::string("test") );

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;
}


StatusCode UseTholdMuon::execute( ) {

  MsgStream log(msgSvc(), name());

  //  std::string fullPath = "/Calib/CAL_TholdMuon/test";
  DataObject *pObject;
  

  m_pCalibDataSvc->retrieveObject(m_path, pObject);

  CalibData::CalTholdMuonCol* pTholdMuon = 0;
  pTholdMuon = dynamic_cast<CalibData::CalTholdMuonCol *> (pObject);
  if (!pTholdMuon) {
    log << MSG::ERROR << "Dynamic cast to CalTholdMuonCol failed" << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNo = pTholdMuon->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO 
        << "Processing new TholdMuon calibration after retrieveObject" 
        << endreq;
    m_ser = newSerNo;
    processNew(pTholdMuon, m_path);
  }
  m_pCalibDataSvc->updateObject(pObject);

  pTholdMuon = 0;
  try {
    pTholdMuon = dynamic_cast<CalibData::CalTholdMuonCol *> (pObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to CalTholdMuonCol after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pTholdMuon->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new TholdMuon after update" 
        << endreq;
    m_ser = newSerNo;
    processNew(pTholdMuon, m_path);
  }

  return StatusCode::SUCCESS;
}


void UseTholdMuon::processNew(CalibData::CalTholdMuonCol* pNew, 
                              const std::string& path) {
  using idents::CalXtalId;
  using CalibData::CalTholdMuon;
  using CalibData::ValSig;
  using CalibData::Xpos;

  bool  done = false;

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Retrieved with path " << path << endreq
      << "Serial #" <<  pNew->getSerNo() << endreq; 
  log << MSG::INFO << "Vstart: " <<  (pNew->validSince()).hour(true)
      << "  Vend: " << (pNew->validTill()).hour(true) << endreq;
  
  if (!done) {
    done = true;
    short iTower = 0;
    short iLayer = 0;
    short iXtal = 1;
    //   No range for this type
    unsigned range = 0;
    unsigned face = 0;
    CalXtalId id(iTower, iLayer, iXtal);
    
    CalibData::RangeBase* pRange = pNew->getRange(id, range, face);
    
    CalTholdMuon* pTholdMuon = dynamic_cast<CalTholdMuon * >(pRange);
    log << MSG::INFO << "For tower = " << iTower << " layer = " << iLayer
        << " xtal = " << iXtal << endreq;

    log << MSG::INFO << "Retrieving TholdMuon values and uncertainties:"  
        << endreq;

    if (!pTholdMuon) {
      log << MSG::INFO << "No calibration data found for this channel" 
          << endreq;
      return;
    }

    const ValSig* pFLE = pTholdMuon->getFLE();
    const ValSig* pFHE = pTholdMuon->getFHE();

    printValSig(&log, pFLE, "FLE: ");
    printValSig(&log, pFHE, "FHE: ");


    const ValSig* pPed = pTholdMuon->getPed(CalXtalId::LEX8);
    if (pPed) {
      printValSig(&log, pPed, "ped for LEX8");
    }
    else {
      log << "No ped for LEX8" << endreq;
    }
    pPed = pTholdMuon->getPed(CalXtalId::LEX1);
    if (pPed) {
      printValSig(&log, pPed, "ped for LEX1");
    }
    else {
      log << "No ped for LEX1" << endreq;
    }

    pPed = pTholdMuon->getPed(CalXtalId::HEX8);
    if (pPed) {
      printValSig(&log, pPed, "ped for HEX8");
    }
    else {
      log << "No ped for HEX8" << endreq;
    }

    pPed = pTholdMuon->getPed(CalXtalId::HEX1);
    if (pPed) {
      printValSig(&log, pPed, "ped for HEX1");
    }
    else {
      log << "No ped for HEX1" << endreq;
    }


  }
}

void UseTholdMuon::printValSigs(MsgStream* log, 
                           const std::vector<CalibData::ValSig>* v, 
                           std::string title) {
  (*log) << MSG::INFO << "ValSig vector " << title << endreq;
  unsigned n = v->size();
  for (unsigned i = 0; i < n; i++) {
    (*log) << MSG::INFO << "val = " << (*v)[i].m_val 
           << "  sig = " << (*v)[i].m_sig << endreq;
  }
  return;
}

void UseTholdMuon::printValSig(MsgStream* log, 
                               const CalibData::ValSig* v, 
                               std::string title) {
  (*log) << MSG::INFO << "ValSig  " << title << endreq;
  (*log) << MSG::INFO << "val = " << v->m_val 
         << "  sig = " << v->m_sig << endreq;
  std::string def = (v->isDefined()) ? "is defined" : "is not defined";
  (*log) << def << endreq;

  return;
}

StatusCode UseTholdMuon::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "Finalize UseTholdMuon "
      << endreq;
  
  return StatusCode::SUCCESS;
}


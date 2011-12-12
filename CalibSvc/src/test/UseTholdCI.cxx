//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibSvc/ICalibPathSvc.h"
#include "CalibData/Cal/CalTholdCI.h"
#include "CalibData/Cal/Xpos.h"
#include "idents/CalXtalId.h"   

/**
   @file UseTholdCI.cxx
   Simple algorithm to test functioning of "the other" TDS, Cal new 
   (Sept. 04) tholdCI calibration data.
*/

  /** 
   @class UseTholdCI

   Algorithm exemplifying retrieval and use of CAL_TholdCI calibration.
*/
class UseTholdCI : public Algorithm {

public:
  UseTholdCI(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper function called by execute
  void processNew(CalibData::CalTholdCICol* pNew, const std::string& path);

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
//static const AlgFactory<UseTholdCI> Factory;
//const IAlgFactory& UseTholdCIFactory = Factory;


UseTholdCI::UseTholdCI(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
  m_ser(-1)
{
  // Declare properties here.

}


StatusCode UseTholdCI::initialize() {
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
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_CAL_TholdCI,
                                  std::string("test") );

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;
}


StatusCode UseTholdCI::execute( ) {

  MsgStream log(msgSvc(), name());

  //  std::string fullPath = "/Calib/CAL_TholdCI/test";
  DataObject *pObject;
  

  m_pCalibDataSvc->retrieveObject(m_path, pObject);

  CalibData::CalTholdCICol* pTholdCI = 0;
  pTholdCI = dynamic_cast<CalibData::CalTholdCICol *> (pObject);
  if (!pTholdCI) {
    log << MSG::ERROR << "Dynamic cast to CalTholdCICol failed" << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNo = pTholdCI->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO 
        << "Processing new TholdCI calibration after retrieveObject" 
        << endreq;
    m_ser = newSerNo;
    processNew(pTholdCI, m_path);
  }
  m_pCalibDataSvc->updateObject(pObject);

  pTholdCI = 0;
  try {
    pTholdCI = dynamic_cast<CalibData::CalTholdCICol *> (pObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to CalTholdCICol after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pTholdCI->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new TholdCI after update" 
        << endreq;
    m_ser = newSerNo;
    processNew(pTholdCI, m_path);
  }

  return StatusCode::SUCCESS;
}


void UseTholdCI::processNew(CalibData::CalTholdCICol* pNew, 
                              const std::string& path) {
  using idents::CalXtalId;
  using CalibData::CalTholdCI;
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
    short iLayer = 2;
    short iXtal = 1;
    // no range for this type
    unsigned range = 0;
    unsigned face = 0;
    CalXtalId id(iTower, iLayer, iXtal);
    
    CalibData::RangeBase* pRange = pNew->getRange(id, range, face);
    
    CalTholdCI* pTholdCI = dynamic_cast<CalTholdCI * >(pRange);
    log << MSG::INFO << "For tower = " << iTower << " layer = " << iLayer
        << " xtal = " << iXtal << endreq;

    log << MSG::INFO << "Retrieving TholdCI values and uncertainties:"  
        << endreq;

    if (!pTholdCI) {
      log << MSG::INFO << "No calibration data found for this channel" 
          << endreq;
      return;
    }

    const ValSig* pFLE = pTholdCI->getFLE();
    const ValSig* pFHE = pTholdCI->getFHE();
    const ValSig* pLAC = pTholdCI->getLAC();

    printValSig(&log, pFLE, "FLE: ");
    printValSig(&log, pFHE, "FHE: ");
    printValSig(&log, pLAC, "LAC: ");


    const ValSig* pPed = pTholdCI->getPed(CalXtalId::LEX8);
    const ValSig* pULD = pTholdCI->getULD(CalXtalId::LEX8);
    if (pPed) {
      printValSig(&log, pPed, "ped for LEX8");
      printValSig(&log, pULD, "ULD for LEX8");
    }
    else {
      log << "No data for LEX8" << endreq;
    }
    pPed = pTholdCI->getPed(CalXtalId::LEX1);
    pULD = pTholdCI->getULD(CalXtalId::LEX1);
    if (pPed) {
      printValSig(&log, pPed, "ped for LEX1");
      printValSig(&log, pULD, "ULD for LEX1");
    }
    else {
      log << "No data for LEX1" << endreq;
    }

    pPed = pTholdCI->getPed(CalXtalId::HEX8);
    pULD = pTholdCI->getULD(CalXtalId::HEX8);
    if (pPed) {
      printValSig(&log, pPed, "ped for HEX8");
      printValSig(&log, pULD, "ULD for HEX8");
    }
    else {
      log << "No data for HEX8" << endreq;
    }

    pPed = pTholdCI->getPed(CalXtalId::HEX1);
    pULD = pTholdCI->getULD(CalXtalId::HEX1);
    if (pPed) {
      printValSig(&log, pPed, "ped for HEX1");
      printValSig(&log, pULD, "ULD for HEX1");
    }
    else {
      log << "No data for HEX1" << endreq;
    }


  }
}

void UseTholdCI::printValSigs(MsgStream* log, 
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

void UseTholdCI::printValSig(MsgStream* log, 
                               const CalibData::ValSig* v, 
                               std::string title) {
  (*log) << MSG::INFO << "ValSig  " << title << endreq;
  (*log) << MSG::INFO << "val = " << v->m_val 
         << "  sig = " << v->m_sig << endreq;
  std::string def = (v->isDefined()) ? "is defined" : "is not defined";
  (*log) << def << endreq;

  return;
}

StatusCode UseTholdCI::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "Finalize UseTholdCI "
      << endreq;
  
  return StatusCode::SUCCESS;
}


//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Acd/AcdCalibGain.h"
#include "CalibData/Acd/AcdCalibPed.h"
#include "idents/AcdId.h"

/**
   @file UseAcd.cxx                         
   Simple algorithm to test functioning of "the other" TDS with
   ACD peds and gains
*/



  /** 
   @class UseAcd

   Algorithm exemplifying retrieval and use of ACD calibrations
*/
class UseAcd : public Algorithm {


public:
  UseAcd(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper functions called by execute
  void processNew(CalibData::AcdCalibGain* pNewGain, 
                  const std::string& pathGain,
                  CalibData::AcdCalibPed* pNewPed, 
                  const std::string& pathPed);

  IDataProviderSvc* m_pCalibDataSvc;
  int               m_serGain;
  int               m_serPed;
};


/// Instantiation of a static factory to create instances of this algorithm
static const AlgFactory<UseAcd> Factory;
const IAlgFactory& UseAcdFactory = Factory;


UseAcd::UseAcd(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
    m_serGain(-1),m_serPed(-1)
{
  // Declare properties here.

}


StatusCode UseAcd::initialize() {
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

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;

}


StatusCode UseAcd::execute( ) {

  MsgStream log(msgSvc(), name());

  std::string fullPathGain = "/Calib/ACD_ElecGain/vanilla";
  std::string fullPathPed = "/Calib/ACD_Ped/vanilla";
  DataObject *pObjectGain;
  DataObject *pObjectPed;
  

  m_pCalibDataSvc->retrieveObject(fullPathGain, pObjectGain);
  m_pCalibDataSvc->retrieveObject(fullPathPed, pObjectPed);

  CalibData::AcdCalibGain* pGains = 0;
  CalibData::AcdCalibPed* pPeds = 0;
  pGains = dynamic_cast<CalibData::AcdCalibGain *> (pObjectGain);
  if (!pGains) {
    log << MSG::ERROR << "Dynamic cast to AcdCalibGain failed" << endreq;
    return StatusCode::FAILURE;
  }

  pPeds = dynamic_cast<CalibData::AcdCalibPed *> (pObjectPed);
  if (!pPeds) {
    log << MSG::ERROR << "Dynamic cast to AcdCalibPed failed" << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNoGain = pGains->getSerNo();
  int newSerNoPed = pPeds->getSerNo();
  if ((newSerNoGain != m_serGain) ||
      (newSerNoPed != m_serPed) )  {
    log << MSG::INFO << "Processing new gains, peds after retrieveObject" 
        << endreq;
    m_serGain = newSerNoGain;
    m_serPed = newSerNoPed;
    processNew(pGains, fullPathGain, pPeds, fullPathPed);
  }
  m_pCalibDataSvc->updateObject(pObjectGain);
  m_pCalibDataSvc->updateObject(pObjectPed);


  pGains = 0;
  pPeds = 0;
  try {
    pGains = dynamic_cast<CalibData::AcdCalibGain *> (pObjectGain);
    pPeds = dynamic_cast<CalibData::AcdCalibPed *> (pObjectPed);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to AcdCalibGain or AcdCalibPed after update failed" 
        << endreq;
    return StatusCode::FAILURE;
  }
  newSerNoGain = pGains->getSerNo();
  newSerNoPed = pPeds->getSerNo();
  if ((newSerNoGain != m_serGain)  ||
      (newSerNoPed != m_serPed)  )   {
    log << MSG::INFO << "Processing new gains, peds after update" 
        << endreq;
    m_serGain = newSerNoGain;
    m_serPed = newSerNoPed;
    processNew(pGains, fullPathGain, pPeds, fullPathPed);
  }

  return StatusCode::SUCCESS;
}


void UseAcd::processNew(CalibData::AcdCalibGain* pNewGain, 
                        const std::string& pathGain,
                        CalibData::AcdCalibPed* pNewPed,
                        const std::string& pathPed) {
  using CalibData::AcdGain;
  using CalibData::AcdPed;
  bool  done = false;

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Retrieved with path " << pathGain << endreq
      << "Serial #" <<  pNewGain->getSerNo() << endreq; 
  log << MSG::INFO << "Vstart: " <<  (pNewGain->validSince()).hours()
      << "  Vend: " << (pNewGain->validTill()).hours() << endreq;
  
  log << MSG::INFO << "And retrieved with path " << pathPed << endreq
      << "Serial #" <<  pNewPed->getSerNo() << endreq; 
  log << MSG::INFO << "Vstart: " <<  (pNewPed->validSince()).hours()
      << "  Vend: " << (pNewPed->validTill()).hours() << endreq;
  
  if (!done) {
    done = true;
    short iFace = 0;
    short iRow = 1;
    short iCol = 2;
    unsigned pmt = 1;
    unsigned range = 0;
    //    CalXtalId id(iTower, iLayer, iXtal);
    // First arg. to constructor is "layer"
    idents::AcdId id(0, iFace, iRow, iCol);
    
    CalibData::RangeBase* pRangeGain = pNewGain->getRange(id, pmt, range);
    CalibData::RangeBase* pRangePed = pNewPed->getRange(id, pmt, range);
    
    AcdGain* pGain = dynamic_cast<AcdGain * >(pRangeGain);
    AcdPed* pPed = dynamic_cast<AcdPed * >(pRangePed);
    log << MSG::INFO << "For face = " << iFace << " row = " << iRow
        << " column = " << iCol << endreq;
    log << MSG::INFO << " pmt = " << pmt 
        << "    range = " << range << endreq;
    
    log << MSG::INFO << " gain = " << pGain->getGain() << endreq;
    log << MSG::INFO << " sigma gain = " << pGain->getSig() << endreq;
    log << MSG::INFO << " ped = " << pPed->getPed() << endreq;
    log << MSG::INFO << " sigma ped = " << pPed->getSig() << endreq;

    iFace = 1;

    idents::AcdId idf1(0, iFace, iRow, iCol);
    
    pRangeGain = pNewGain->getRange(idf1, pmt, range);
    pRangePed = pNewPed->getRange(idf1, pmt, range);
    
    pGain = dynamic_cast<AcdGain * >(pRangeGain);
    pPed = dynamic_cast<AcdPed * >(pRangePed);
    log << MSG::INFO << "For face = " << iFace << " row = " << iRow
        << " column = " << iCol << endreq;
    log << MSG::INFO << " pmt = " << pmt 
        << "    range = " << range << endreq;
    
    log << MSG::INFO << " gain = " << pGain->getGain() << endreq;
    log << MSG::INFO << " sigma gain = " << pGain->getSig() << endreq;
    log << MSG::INFO << " ped = " << pPed->getPed() << endreq;
    log << MSG::INFO << " sigma ped = " << pPed->getSig() << endreq;


  }
}

StatusCode UseAcd::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "          Finalize UseAcd "
      << endreq;
  
  return StatusCode::SUCCESS;
}


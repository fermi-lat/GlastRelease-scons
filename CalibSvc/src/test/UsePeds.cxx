//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Cal/CalCalibPed.h"
#include "CalibData/CalibTime.h"
#include "idents/CalXtalId.h"                // shouldn't be necessary

/**
   @file UsePeds.cxx
   Simple algorithm to test functioning of "the other" TDS, Cal pedestals data
*/



  /** 
   @class UsePeds

   Algorithm exemplifying retrieval and use of Calorimeter pedestal calibration
*/
class UsePeds : public Algorithm {


public:
  UsePeds(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper function called by execute
  void processNew(CalibData::CalCalibPed* pNew, const std::string& path);

  IDataProviderSvc* m_pCalibDataSvc;
  int               m_ser;
};


/// Instantiation of a static factory to create instances of this algorithm
static const AlgFactory<UsePeds> Factory;
const IAlgFactory& UsePedsFactory = Factory;


UsePeds::UsePeds(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
  m_ser(-1)
{
  // Declare properties here.

}


StatusCode UsePeds::initialize() {
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


StatusCode UsePeds::execute( ) {

  MsgStream log(msgSvc(), name());

  //  SmartDataPtr<CalibData::CalibTest1> test1(m_pCalibDataSvc,
  //                                        CalibData::Test_Gen);
  //  CalibData::CalibTest1* test1 = 
  //    SmartDataPtr<CalibData::CalibTest1>(m_pCalibDataSvc, CalibData::Test_Gen);
  
  std::string fullPath = "/Calib/CAL_Ped/ideal";
  DataObject *pObject;
  

  m_pCalibDataSvc->retrieveObject(fullPath, pObject);

  CalibData::CalCalibPed* pPeds = 0;
  pPeds = dynamic_cast<CalibData::CalCalibPed *> (pObject);
  if (!pPeds) {
    log << MSG::ERROR << "Dynamic cast to CalCalibPed failed" << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNo = pPeds->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new pedestals after retrieveObject" 
        << endreq;
    m_ser = newSerNo;
    processNew(pPeds, fullPath);
  }
  m_pCalibDataSvc->updateObject(pObject);

  pPeds = 0;
  try {
    pPeds = dynamic_cast<CalibData::CalCalibPed *> (pObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to CalCalibPed after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pPeds->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new pedestals after update" 
        << endreq;
    m_ser = newSerNo;
    processNew(pPeds, fullPath);
  }

  return StatusCode::SUCCESS;
}


void UsePeds::processNew(CalibData::CalCalibPed* pNew, 
                              const std::string& path) {
  using idents::CalXtalId;
  using CalibData::Ped;
  bool  done = false;

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Retrieved with path " << path << endreq
      << "Serial #" <<  pNew->getSerNo() << endreq; 
  log << MSG::INFO << "Vstart: " <<  (pNew->validSince()).hours()
      << "  Vend: " << (pNew->validTill()).hours() << endreq;

  std::string vStart = pNew->getValidStart()->getString();
  std::string vEnd = pNew->getValidEnd()->getString();
  log << MSG::INFO << "(Vstart, Vend) as ascii from CalibTime objects: " 
      << endreq;
  log << MSG::INFO << "(" << vStart << ", " << vEnd << ")" << endreq;

  
  if (!done) {
    done = true;
    short iTower = 0;
    short iLayer = 0;
    short iXtal = 2;
    unsigned range = 2;
    unsigned face = 0;
    CalXtalId id(iTower, iLayer, iXtal);
    
    CalibData::RangeBase* pRange = pNew->getRange(id, range, face);
    
    Ped* pPed = dynamic_cast<Ped * >(pRange);
    log << MSG::INFO << "For tower = " << iTower << " layer = " << iLayer
        << " xtal = " << iXtal << endreq;
    log << MSG::INFO << "    range = " << range 
        << " face = " << face << endreq;
    
    log << MSG::INFO << "Averaged ped = " << pPed->getAvr() << endreq;
    log << MSG::INFO << "       sigma = " << pPed->getSig() << endreq;
    log << MSG::INFO << "  cos angle = " << pPed->getCosAngle() << endreq;

    /*      Try another tower */
    iTower++;
    id = CalXtalId(iTower, iLayer, iXtal);
    
    pRange = pNew->getRange(id, range, face);
    
    pPed = dynamic_cast<Ped * >(pRange);
    log << MSG::INFO << "For tower = " << iTower << " layer = " << iLayer
        << " xtal = " << iXtal << endreq;
    log << MSG::INFO << "    range = " << range 
        << " face = " << face << endreq;
    
    log << MSG::INFO << "Averaged ped = " << pPed->getAvr() << endreq;
    log << MSG::INFO << "       sigma = " << pPed->getSig() << endreq;
    log << MSG::INFO << "  cos angle = " << pPed->getCosAngle() << endreq;

  }
}

StatusCode UsePeds::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "          Finalize UsePeds "
      << endreq;
  
  return StatusCode::SUCCESS;
}


//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Cal/CalCalibGain.h"
#include "idents/CalXtalId.h"                // shouldn't be necessary

/**
   @file UseGains.cxx                         
   Simple algorithm to test functioning of "the other" TDS, Cal gains data
*/



  /** 
   @class UseGains

   Algorithm exemplifying retrieval and use of Calorimeter gain calibration
*/
class UseGains : public Algorithm {


public:
  UseGains(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper function called by execute
  void processNew(CalibData::CalCalibGain* pNew, const std::string& path);

  IDataProviderSvc* m_pCalibDataSvc;
  int               m_ser;
};


/// Instantiation of a static factory to create instances of this algorithm
static const AlgFactory<UseGains> Factory;
const IAlgFactory& UseGainsFactory = Factory;


UseGains::UseGains(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
  m_ser(-1)
{
  // Declare properties here.

}


StatusCode UseGains::initialize() {
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


StatusCode UseGains::execute( ) {

  MsgStream log(msgSvc(), name());

  //  SmartDataPtr<CalibData::CalibTest1> test1(m_pCalibDataSvc,
  //                                        CalibData::Test_Gen);
  //  CalibData::CalibTest1* test1 = 
  //    SmartDataPtr<CalibData::CalibTest1>(m_pCalibDataSvc, CalibData::Test_Gen);
  
  std::string fullPath = "/Calib/CAL_ElecGain/vanilla";
  DataObject *pObject;
  

  m_pCalibDataSvc->retrieveObject(fullPath, pObject);

  CalibData::CalCalibGain* pGains = 0;
  pGains = dynamic_cast<CalibData::CalCalibGain *> (pObject);
  if (!pGains) {
    log << MSG::ERROR << "Dynamic cast to CalCalibGain failed" << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNo = pGains->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new gains after retrieveObject" 
        << endreq;
    m_ser = newSerNo;
    processNew(pGains, fullPath);
  }
  m_pCalibDataSvc->updateObject(pObject);

  pGains = 0;
  try {
    pGains = dynamic_cast<CalibData::CalCalibGain *> (pObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to CalCalibGain after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pGains->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new gains after update" 
        << endreq;
    m_ser = newSerNo;
    processNew(pGains, fullPath);
  }

  return StatusCode::SUCCESS;
}


void UseGains::processNew(CalibData::CalCalibGain* pNew, 
                              const std::string& path) {
  using idents::CalXtalId;
  using CalibData::Gain;
  bool  done = false;

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Retrieved with path " << path << endreq
      << "Serial #" <<  pNew->getSerNo() << endreq; 
  log << MSG::INFO << "Vstart: " <<  (pNew->validSince()).hours()
      << "  Vend: " << (pNew->validTill()).hours() << endreq;
  
  if (!done) {
    done = true;
    short iTower = 0;
    short iLayer = 0;
    short iXtal = 2;
    unsigned range = 2;
    unsigned face = 0;
    CalXtalId id(iTower, iLayer, iXtal);
    
    CalibData::RangeBase* pRange = pNew->getRange(id, range, face);
    
    Gain* pGain = dynamic_cast<Gain * >(pRange);
    log << MSG::INFO << "For tower = " << iTower << " layer = " << iLayer
        << " xtal = " << iXtal << endreq;
    log << MSG::INFO << "    range = " << range 
        << " face = " << face << endreq;
    
    log << MSG::INFO << "Averaged gain = " << pGain->getGain() << endreq;
    log << MSG::INFO << "       sigma = " << pGain->getSig() << endreq;
  }
}

StatusCode UseGains::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "          Finalize UseGains "
      << endreq;
  
  return StatusCode::SUCCESS;
}


//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibSvc/ICalibRootSvc.h"
#include "CalibSvc/ICalibPathSvc.h"
#include "CalibData/Cal/CalCalibGain.h"
#include "idents/CalXtalId.h"                // shouldn't be necessary
#include "facilities/commonUtilities.h"

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
  ICalibPathSvc*    m_pCalibPathSvc;
  ICalibRootSvc*    m_pRootSvc;
  int               m_ser;
  int               m_serR;
  std::string       m_path;
  std::string       m_pathR;
  bool              m_didWrite;
  std::string       m_outName;
  bool              m_doWrite;
};


/// Instantiation of a static factory to create instances of this algorithm
static const AlgFactory<UseGains> Factory;
const IAlgFactory& UseGainsFactory = Factory;


UseGains::UseGains(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
    m_ser(-1), m_serR(-1), m_didWrite(false)
{
  // Declare properties here.
  declareProperty("OutName", m_outName="CalGains.root");
  declareProperty("DoWrite", m_doWrite= false);

}


StatusCode UseGains::initialize() {
  StatusCode sc;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initialize()" << endreq;


  sc = service("CalibDataSvc", m_pCalibDataSvc, true);

  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Could not get IDataProviderSvc interface of CalibDataSvc" 
	<< endreq;
    return sc;
  }

#ifdef SCons
  // define CALIBUTILROOT.  It's embedded in metadata for our calib
  std::string calibUtilRoot = 
    facilities::commonUtilities::getPackagePath("calibUtil");
  facilities::commonUtilities::setEnvironment("CALIBUTILROOT", calibUtilRoot);
#endif
  sc = service("CalibDataSvc", m_pCalibPathSvc, true);

  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Could not get ICalibPathSvc interface of CalibDataSvc" 
	<< endreq;
    return sc;
  }

  m_path = 
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_CAL_ElecGain,
                                  std::string("ideal") );

  m_pathR = 
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_CAL_ElecGain,
                                  std::string("RootBeer") );

  sc = service("CalibRootCnvSvc", m_pRootSvc, true);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Could not get ICalibRootSvc interface of CalibRootCnvSvc" 
	<< endreq;
    return sc;
  }

  // Get properties from the JobOptionsSvc
  sc = setProperties();

  // If write isn't requested, just pretend we already did it.
  m_didWrite = !m_doWrite;

  // Probably should be set by job options.
  //  m_outName = "CalibGains.root";

  return StatusCode::SUCCESS;

}


StatusCode UseGains::execute( ) {

  MsgStream log(msgSvc(), name());

  //  std::string fullPath = "/Calib/CAL_ElecGain/ideal";
  DataObject *pObject;

  m_pCalibDataSvc->retrieveObject(m_path, pObject);

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
    processNew(pGains, m_path);
    if (!m_didWrite) {
      StatusCode writeStatus = 
        m_pRootSvc->writeToRoot(m_outName, m_path);
      log << "Status returned from writing " << m_outName << " is " 
          << writeStatus << endreq;
      m_didWrite = true;
    }
      
  }
  // Same thing with RootBeer
  //  std::string fullPathR = "/Calib/CAL_ElecGain/RootBeer";
  DataObject *pObjectR;
  m_pCalibDataSvc->retrieveObject(m_pathR, pObjectR);

  CalibData::CalCalibGain* pGainsR = 0;
  pGainsR = dynamic_cast<CalibData::CalCalibGain *> (pObjectR);
  if (!pGainsR) {
    log << MSG::ERROR << "Dynamic cast to CalCalibGain (RootBeer) failed" 
        << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNoR = pGainsR->getSerNo();
  if (newSerNoR != m_serR) {
    log << MSG::INFO << "Processing new RootBeer gains after retrieveObject" 
        << endreq;
    m_serR = newSerNoR;
    processNew(pGainsR, m_pathR);
  }


  //////


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
    processNew(pGains, m_path);
    if (!m_didWrite) {
      StatusCode writeStatus = 
        m_pRootSvc->writeToRoot(m_outName, m_path);
      log << "Status returned from writing " << m_outName << " is " 
          << writeStatus << endreq;
      m_didWrite = true;
    }
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


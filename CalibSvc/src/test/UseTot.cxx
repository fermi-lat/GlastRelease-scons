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
#include "CalibData/Tkr/TkrTot.h"

/**
   @file UseTot.cxx                         
   Simple algorithm to test functioning of "the other" TDS, TKR 
   charge injection ToT data
*/

  /** 
   @class UseTot

   Algorithm exemplifying retrieval and use of Tracker ToT calibration
*/
class UseTot : public Algorithm {


public:
  UseTot(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper function called by execute
  void processNew(CalibData::TkrTotCol* pNew, const std::string& path);

  IDataProviderSvc* m_pCalibDataSvc;
  ICalibRootSvc*    m_pRootSvc;
  ICalibPathSvc*    m_pCalibPathSvc;
  std::string       m_path;
  int               m_ser;
  int               m_serR;
};


/// Instantiation of a static factory to create instances of this algorithm
//static const AlgFactory<UseTot> Factory;
//const IAlgFactory& UseTotFactory = Factory;
DECLARE_ALGORITHM_FACTORY(UseTot);


UseTot::UseTot(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
    m_ser(-1), m_serR(-1)
{
  // Declare properties here.

}


StatusCode UseTot::initialize() {
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

  sc = service("CalibRootCnvSvc", m_pRootSvc, true);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Could not get ICalibRootSvc interface of CalibRootCnvSvc" 
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
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_TKR_TOTSignal,
                                  std::string("RootBeer") );

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;

}


StatusCode UseTot::execute( ) {

  using CalibData::TkrTotCol;
  MsgStream log(msgSvc(), name());

  //  std::string fullPath = "/Calib/TKR_TOTSignal/RootBeer";
  DataObject *pObject;

  m_pCalibDataSvc->retrieveObject(m_path, pObject);

  TkrTotCol* pTotCol = 0;
  pTotCol = dynamic_cast<TkrTotCol *> (pObject);
  if (!pTotCol) {
    log << MSG::ERROR << "Dynamic cast to TkrTotCol (RootBeer) failed" 
        << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNo = pTotCol->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new TKR_TOTSignal after retrieveObject" 
        << endreq;
    m_ser = newSerNo;
    processNew(pTotCol, m_path);
      
  }


  m_pCalibDataSvc->updateObject(pObject);

  pTotCol = 0;
  try {
    pTotCol = dynamic_cast<TkrTotCol *> (pObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to TkrTotCol after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pTotCol->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new ToT collection after update" 
        << endreq;
    m_ser = newSerNo;
    processNew(pTotCol, m_path);
  }

  return StatusCode::SUCCESS;
}


void UseTot::processNew(CalibData::TkrTotCol* pNew, 
                              const std::string& path) {
  using idents::TkrId;
  using CalibData::TkrTotCol;
  using CalibData::TkrTotUni;
  bool  done = false;

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Retrieved with path " << path << endreq
      << "Serial #" <<  pNew->getSerNo() << endreq; 
  log << MSG::INFO << "Vstart: " <<  (pNew->validSince()).hour(true)
      << "  Vend: " << (pNew->validTill()).hour(true) << endreq;
  
  unsigned tray = 0;
  bool top = true;
  while (!done) {
    unsigned towerX = 0;
    unsigned towerY = 0;
    //    unsigned tray = 3;
    //    bool top = true;
    TkrId id(towerX, towerY, tray, top);
        if ((tray == 0) && (top)) {
      const std::string* hwserial = pNew->getHwserial(towerY, towerX);

      log << MSG::INFO << "hw serial for tower x = " << towerX
          << " and tower y = " << towerY << " is " << *hwserial << std::endl;
    }

    unsigned iStrip = 27;

    const CalibData::TkrTotStrip* pInfo = pNew->getStripInfo(id, iStrip);
    
    log << MSG::INFO << "For tower X,Y = (" << towerX << ", " << towerY
        << ") tray = " << tray << " top = " << top << " strip = "
        << iStrip << endreq;
    if (!pInfo) {
      log << MSG::INFO << "No info found " << endreq;
    }
    else {
      log << MSG::INFO << "Strip id  = " << pInfo->getStripId() << endreq;
      log << MSG::INFO << "Slope     = " << pInfo->getSlope() << endreq;
      log << MSG::INFO << "Intercept = " << pInfo->getIntercept() << endreq;
      log << MSG::INFO << "Quad      = " << pInfo->getQuad() << endreq;
      log << MSG::INFO << "Chi2      = " << pInfo->getChi2() << endreq;
      log << MSG::INFO << "Deg. fr.  = " << pInfo->getDf() << endreq;
    }
    if (top) tray++;
    else if (tray == 18) done = true;
    top = !top;
  }
}

StatusCode UseTot::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "          Finalize UseTot "
      << endreq;
  
  return StatusCode::SUCCESS;
}


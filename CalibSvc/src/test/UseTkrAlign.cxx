//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibSvc/ICalibPathSvc.h"
#include "CalibData/Tkr/TkrTowerAlignCalib.h"
#include "CalibData/Tkr/TkrInternalAlignCalib.h"
#include "CalibData/CalibTime.h"
#include "idents/TkrId.h"

/**
   @file UseTkr.cxx
   Simple algorithm to test functioning of "the other" TDS, 
   Tkr alignment (both inter-tower and internal to towers)
*/



  /** 
   @class UseTkrAlign

   Algorithm exemplifying retrieval and use of Tkr calibration quantities
*/

namespace {
  void outputAlign(MsgStream& log, const CLHEP::Hep3Vector& disp, 
                   const CLHEP::Hep3Vector& rot) {
    double x = disp.x(), y = disp.y(), z = disp.z();
    log << "  disp=(" << x << ", " << y << ", " << z << ") " << endreq;
    x =  rot.x(); y = rot.y(); z = rot.z();
    log << "  rot=(" << x << ", " << y << ", " << z  << ") ";
    log << endreq << endreq;
  }
}


class UseTkrAlign : public Algorithm {


public:
  UseTkrAlign(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper functions called by execute
  void processInter(CalibData::TkrTowerAlignCalib* pNew, 
                    const std::string& path);
  void processIntra(CalibData::TkrInternalAlignCalib* pNew, 
                    const std::string& path);

  IDataProviderSvc* m_pCalibDataSvc;
  ICalibPathSvc*    m_pCalibPathSvc;
  std::string       m_pathInter;
  std::string       m_pathIntra;
  int               m_towerSer;
  int               m_internalSer;
};


/// Instantiation of a static factory to create instances of this algorithm
//static const AlgFactory<UseTkrAlign> Factory;
//const IAlgFactory& UseTkrAlignFactory = Factory;
DECLARE_ALGORITHM_FACTORY(UseTkrAlign);

UseTkrAlign::UseTkrAlign(const std::string&  name, 
               ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
    m_towerSer(-1), m_internalSer(-1)
{
  // Declare properties here.

}


StatusCode UseTkrAlign::initialize() {
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

  m_pathInter = 
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_TKR_TowerAlign,
                                  std::string("vanilla") );
  m_pathIntra = 
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_TKR_InternalAlign,
                                  std::string("vanilla") );

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;

}


StatusCode UseTkrAlign::execute( ) {

  MsgStream log(msgSvc(), name());


  DataObject *pObject;
  

  m_pCalibDataSvc->retrieveObject(m_pathInter, pObject);

  CalibData::TkrTowerAlignCalib* pTowerAlign = 0;
  pTowerAlign = dynamic_cast<CalibData::TkrTowerAlignCalib *> (pObject);
  if (!pTowerAlign) {
    log << MSG::ERROR << "Dynamic cast to TkrTowerAlignCalib failed" << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNo = pTowerAlign->getSerNo();
  if (newSerNo != m_towerSer) {
    log << MSG::INFO << "Processing new tower align after retrieveObject" 
        << endreq;
    m_towerSer = newSerNo;
    processInter(pTowerAlign, m_pathInter);
  }
  m_pCalibDataSvc->updateObject(pObject);

  pTowerAlign = 0;
  try {
    pTowerAlign = dynamic_cast<CalibData::TkrTowerAlignCalib *> (pObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to TkrTowerAlignCalib after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pTowerAlign->getSerNo();
  if (newSerNo != m_towerSer) {
    log << MSG::INFO << "Processing new tower align after update" 
        << endreq;
    m_towerSer = newSerNo;
    processInter(pTowerAlign, m_pathInter);
  }


  DataObject *pObject2;
  CalibData::TkrInternalAlignCalib* pInternalAlign = 0;
  m_pCalibDataSvc->retrieveObject(m_pathIntra, pObject2);
  pInternalAlign = dynamic_cast<CalibData::TkrInternalAlignCalib *> (pObject2);
  if (!pInternalAlign) {
    log << MSG::ERROR << "Dynamic cast to TkrInternalAlignCalib failed" 
        << endreq;
    return StatusCode::FAILURE;
  }

  newSerNo = pInternalAlign->getSerNo();
  if (newSerNo != m_internalSer) {
    log << MSG::INFO << "Processing new internal align after retrieveObject" 
        << endreq;
    m_internalSer = newSerNo;
    processIntra(pInternalAlign, m_pathIntra);
  }
  m_pCalibDataSvc->updateObject(pObject2);

  pInternalAlign = 0;
  try {
    pInternalAlign = dynamic_cast<CalibData::TkrInternalAlignCalib *> (pObject2);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to TkrInternalAlignCalib after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pInternalAlign->getSerNo();
  if (newSerNo != m_internalSer) {
    log << MSG::INFO << "Processing new internal align after update" 
        << endreq;
    m_internalSer = newSerNo;
    processIntra(pInternalAlign, m_pathIntra);
  }

  return StatusCode::SUCCESS;
}


void UseTkrAlign::processInter(CalibData::TkrTowerAlignCalib* pNew, 
                              const std::string& path) {
  bool  done = false;

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Retrieved with path " << path << endreq
      << "Serial #" <<  pNew->getSerNo() << endreq; 
  log << MSG::INFO << "Vstart: " <<  (pNew->validSince()).hour(true)
      << "  Vend: " << (pNew->validTill()).hour(true) << endreq;

  std::string vStart = pNew->getValidStart()->getString();
  std::string vEnd = pNew->getValidEnd()->getString();
  log << MSG::INFO << "(Vstart, Vend) as ascii from CalibTime objects: " 
      << endreq;
  log << MSG::INFO << "(" << vStart << ", " << vEnd << ")" << endreq;

  CLHEP::Hep3Vector disp, rot;

  if (!done) {
    done = true;
    for (unsigned ix = 0; ix < 16; ix++) {
      pNew->getTowerAlign(ix, disp, rot);
      log << MSG::INFO << "Tower alignment for tower #" << ix << ":" 
          << endreq;
      outputAlign(log, disp, rot);
    }
  }
}

void UseTkrAlign::processIntra(CalibData::TkrInternalAlignCalib* pNew, 
                              const std::string& path) {
  bool  done = false;

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Retrieved with path " << path << endreq
      << "Serial #" <<  pNew->getSerNo() << endreq; 
  log << MSG::INFO << "Vstart: " <<  (pNew->validSince()).hour(true)
      << "  Vend: " << (pNew->validTill()).hour(true) << endreq;

  std::string vStart = pNew->getValidStart()->getString();
  std::string vEnd = pNew->getValidEnd()->getString();
  log << MSG::INFO << "(Vstart, Vend) as ascii from CalibTime objects: " 
      << endreq;
  log << MSG::INFO << "(" << vStart << ", " << vEnd << ")" << endreq;

  
  if (!done) {
    // fetch some random values
    done = true;
    unsigned tower = 1;
    unsigned tray = 1;
    CLHEP::Hep3Vector disp, rot;
    StatusCode ok = pNew->getTrayAlign(tower, tray, disp, rot);
    log << "Constants for tower" << tower << ", tray " << tray << ": " 
        << endreq;
    outputAlign(log, disp, rot);

    tray = 3;
    ok = pNew->getTrayAlign(tower, tray, disp, rot);
    log << "Constants for tower" << tower << ", tray " << tray << ": " 
        << endreq;
    outputAlign(log, disp, rot);
    
    tray = 5;
    unsigned face = 0;
    ok = pNew->getFaceAlign(tower, tray, face, disp, rot);
    log << "Constants for tower" << tower << ", tray " << tray 
        << "face " << face << ": " 
        << endreq;
    outputAlign(log, disp, rot);
    
    tower = 6; tray=16; face=0;
    unsigned ladder = 1;
    unsigned wafer = 2;

    ok = pNew->getLadderAlign(tower, tray, face, ladder, disp, rot);
    log << "Constants for tower" << tower << ", tray " << tray 
        << "face " << face << "ladder " << ladder << ": " 
        << endreq;
    outputAlign(log, disp, rot);

    ok = pNew->getWaferAlign(tower, tray, face, ladder, wafer, disp, rot);
    log << "Constants for tower" << tower << ", tray " << tray 
        << "face " << face << "ladder " << ladder << "wafer " << wafer 
        << ": " 
        << endreq;
    outputAlign(log, disp, rot);

  }
}



StatusCode UseTkrAlign::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "          Finalize UseTkrAlign "
      << endreq;
  
  return StatusCode::SUCCESS;
}


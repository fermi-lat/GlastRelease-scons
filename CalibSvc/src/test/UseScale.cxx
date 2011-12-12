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
#include "CalibData/Tkr/TkrScale.h"

/**
   @file UseScale.cxx                         
   Simple algorithm to test functioning of "the other" TDS, TKR 
   charge injection Scale data
*/

  /** 
   @class UseScale

   Algorithm exemplifying retrieval and use of Tracker charge scale calibration
*/
class UseScale : public Algorithm {


public:
  UseScale(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper function called by execute
  void processNew(CalibData::TkrScaleCol* pNew, const std::string& path);

  void outputStrip(unsigned towerX, unsigned towerY, unsigned tray, bool top,
                           CalibData::TkrScaleObj& info, MsgStream& log);

  IDataProviderSvc* m_pCalibDataSvc;
  ICalibRootSvc*    m_pRootSvc;
  ICalibPathSvc*    m_pCalibPathSvc;
  std::string       m_path;
  std::string       m_pathR;  // unused currently
  int               m_ser;
  int               m_serR;   // unused currently
};


/// Instantiation of a static factory to create instances of this algorithm
//static const AlgFactory<UseScale> Factory;
//const IAlgFactory& UseScaleFactory = Factory;
DECLARE_ALGORITHM_FACTORY(UseScale);


UseScale::UseScale(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
    m_ser(-1), m_serR(-1)
{
  // Declare properties here.

}


StatusCode UseScale::initialize() {
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
	<< "Could not get IDataProviderSvc interface of CalibDataSvc" 
	<< endreq;
    return sc;
  }

  m_path = 
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_TKR_ChargeScale,
                                  std::string("RootBeer") );

  // Get properties from the JobOptionsSvc
  sc = setProperties();

  return StatusCode::SUCCESS;

}


StatusCode UseScale::execute( ) {

  using CalibData::TkrScaleCol;
  MsgStream log(msgSvc(), name());

  //  std::string fullPath = "/Calib/TKR_ChargeScale/RootBeer";
  DataObject *pObject;

  m_pCalibDataSvc->retrieveObject(m_path, pObject);

  TkrScaleCol* pScaleCol = 0;
  pScaleCol = dynamic_cast<TkrScaleCol *> (pObject);
  if (!pScaleCol) {
    log << MSG::ERROR << "Dynamic cast to TkrScaleCol (RootBeer) failed" 
        << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNo = pScaleCol->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new TKR_ChargeScale after retrieveObject" 
        << endreq;
    m_ser = newSerNo;
    processNew(pScaleCol, m_path);
      
  }


  m_pCalibDataSvc->updateObject(pObject);

  pScaleCol = 0;
  try {
    pScaleCol = dynamic_cast<TkrScaleCol *> (pObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to TkrScaleCol after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pScaleCol->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new Scale collection after update" 
        << endreq;
    m_ser = newSerNo;
    processNew(pScaleCol, m_path);
  }

  return StatusCode::SUCCESS;
}


void UseScale::processNew(CalibData::TkrScaleCol* pNew, 
                              const std::string& path) {
  using idents::TkrId;
  using CalibData::TkrScaleCol;
  using CalibData::TkrScaleUni;
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
    TkrId id(towerX, towerY, tray, top);
    if ((tray == 0) && (top)) {
      const std::string* hwserial = pNew->getHwserial(towerY, towerX);

      log << MSG::INFO << "hw serial for tower x = " << towerX
          << " and tower y = " << towerY << " is " << *hwserial << std::endl;
    }
    for (unsigned iStrip = 5; iStrip < 1000; iStrip += 65) {

      CalibData::TkrScaleObj info(pNew->getStripInfo(id, iStrip));
      outputStrip(towerX, towerY, tray, top, info, log);

    }
    if (top) tray++;
    else if (tray == 18) done = true;
    top = !top;
  }
}

void UseScale::outputStrip(unsigned towerX, unsigned towerY, unsigned tray,
                           bool top,
                           CalibData::TkrScaleObj& info, MsgStream& log) {
  unsigned iStrip = info.getId();
  log << MSG::INFO << "For tower X,Y = (" << towerX << ", " << towerY
      << ") tray = " << tray << " top = " << top << " strip = "
      << iStrip << endreq;
  if (info.getId() < 0) {
    log << MSG::INFO << "No info found " << endreq;
  }
  else {
    log << MSG::INFO << "Strip id  = " << info.getId() << endreq;
    log << MSG::INFO << "Scale     = " << info.getScale() << endreq;
    log << MSG::INFO << "Error     = " << info.getError() << endreq;
    log << MSG::INFO << "Chi2      = " << info.getChi2() << endreq;
    log << MSG::INFO << "Deg. fr.  = " << info.getDf() << endreq;
  }
}
StatusCode UseScale::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "          Finalize UseScale "
      << endreq;
  
  return StatusCode::SUCCESS;
}


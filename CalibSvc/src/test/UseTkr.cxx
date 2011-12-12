//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibSvc/ICalibPathSvc.h"
#include "CalibData/Tkr/TkrSplitsCalib.h"
#include "CalibData/CalibTime.h"
#include "idents/TkrId.h"

/**
   @file UseTkr.cxx
   Simple algorithm to test functioning of "the other" TDS, 
   Tkr splits data
*/



  /** 
   @class UseTkr

   Algorithm exemplifying retrieval and use of Tkr calibration quantities
*/
class UseTkr : public Algorithm {


public:
  UseTkr(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper function called by execute
  void processNew(CalibData::TkrSplitsCalib* pNew, const std::string& path);

  IDataProviderSvc* m_pCalibDataSvc;
  ICalibPathSvc*    m_pCalibPathSvc;
  std::string       m_path;
  int               m_ser;
};


/// Instantiation of a static factory to create instances of this algorithm
//static const AlgFactory<UseTkr> Factory;
//const IAlgFactory& UseTkrFactory = Factory;
DECLARE_ALGORITHM_FACTORY(UseTkr);


UseTkr::UseTkr(const std::string&  name, 
               ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
  m_ser(-1)
{
  // Declare properties here.

}


StatusCode UseTkr::initialize() {
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
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_TKR_Splits,
                                  std::string("test") );

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;

}


StatusCode UseTkr::execute( ) {

  MsgStream log(msgSvc(), name());


  // std::string fullPath = "/Calib/TKR_Splits/test";
  DataObject *pObject;
  

  m_pCalibDataSvc->retrieveObject(m_path, pObject);

  CalibData::TkrSplitsCalib* pSplits = 0;
  pSplits = dynamic_cast<CalibData::TkrSplitsCalib *> (pObject);
  if (!pSplits) {
    log << MSG::ERROR << "Dynamic cast to TkrSplitsCalib failed" << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNo = pSplits->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new splitss after retrieveObject" 
        << endreq;
    m_ser = newSerNo;
    processNew(pSplits, m_path);
  }
  m_pCalibDataSvc->updateObject(pObject);

  pSplits = 0;
  try {
    pSplits = dynamic_cast<CalibData::TkrSplitsCalib *> (pObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to TkrSplitsCalib after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pSplits->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new splits after update" 
        << endreq;
    m_ser = newSerNo;
    processNew(pSplits, m_path);
  }

  return StatusCode::SUCCESS;
}


void UseTkr::processNew(CalibData::TkrSplitsCalib* pNew, 
                              const std::string& path) {
  using idents::TkrId;
  using CalibData::TkrSplit;
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
    done = true;
    unsigned iRow = 0;
    unsigned iCol = 0;
    unsigned iTray = 1;


    TkrId idTop(iCol, iRow, iTray, true);
    TkrId idBot(iCol, iRow, iTray, false);
    
    CalibData::RangeBase* pTop = pNew->getChannel(idTop);
    CalibData::RangeBase* pBot = pNew->getChannel(idBot);
    
    TkrSplit* pSplitTop = dynamic_cast<TkrSplit * >(pTop);
    TkrSplit* pSplitBot = dynamic_cast<TkrSplit * >(pBot);

    log << MSG::INFO << "For tower row,col = (" << iRow 
        << ", " << iCol << ") tray " << iTray << "top splits are:  "
        << pSplitTop->getLow() << ", " << pSplitTop->getHigh() << endreq;

    log << MSG::INFO << "For tower row,col = (" << iRow 
        << ", " << iCol << ") tray " << iTray << "bot splits are:  "
        << pSplitBot->getLow() << ", " << pSplitBot->getHigh() << endreq;


    /*      Try another tower  -- someday. This file has only one*/

  }
}

StatusCode UseTkr::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "          Finalize UseTkr "
      << endreq;
  
  return StatusCode::SUCCESS;
}


//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Tkr/BadStrips.h"


/**
   @file UseBadStrips.cxx
   Simple algorithm to test functioning of "the other" TDS, bad strips data

   Two classes are defined here.  UseBadStrips is an algorithm which
   accesses hot and dead tracker strip calibration data in what is
   expected to be a standard manner.  The second class, BadVisitor,
   is needed by UseBadStrips to explore the data since, for most
   quantities in the data, the visitor interface is the only way
   to get at them.  
*/



/**
     @class BadVisitor

   Minimal class derived from CalibData::BadStripsVisitor to
   check out BadStrips visitor interface.
*/
class BadVisitor : public CalibData::BadStripsVisitor {
public:
  BadVisitor(MsgStream* log=0) : m_log(log) {}

  void setLog(MsgStream* log) {m_log = log;}

  virtual CalibData::eVisitorRet badTower(unsigned int row, unsigned int col,
                                          int badness);

  virtual CalibData::eVisitorRet badPlane(unsigned int row, unsigned int col, 
                                 unsigned int tray, bool top,
                                 int badness, 
                                 const CalibData::StripCol& strips);
private:
  MsgStream* m_log;
  // Normally would have more here, but for testing purposes all
  // we do is print out stuff from BadStrips as we encounter it
};

/** 
   @class UseBadStrips

   Algorithm exemplifying retrieval and use of calibration bad strips
*/
class UseBadStrips : public Algorithm {


public:
  UseBadStrips(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper function called by execute
  void processNew(CalibData::BadStrips* pNew, const std::string& path);

  IDataProviderSvc* m_pCalibDataSvc;
  int               m_serHot;
  int               m_serDead;
  BadVisitor*       m_visitor;
};


/// Instantiation of a static factory to create instances of this algorithm
static const AlgFactory<UseBadStrips> Factory;
const IAlgFactory& UseBadStripsFactory = Factory;


UseBadStrips::UseBadStrips(const std::string&  name, 
		    ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
  m_serHot(-1), m_serDead(-1)
{
  // Declare properties here.
  m_visitor = new BadVisitor;
}


StatusCode UseBadStrips::initialize() {
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


StatusCode UseBadStrips::execute( ) {

  MsgStream log(msgSvc(), name());
  m_visitor->setLog(&log);

  //  SmartDataPtr<CalibData::CalibTest1> test1(m_pCalibDataSvc,
  //                                        CalibData::Test_Gen);
  //  CalibData::CalibTest1* test1 = 
  //    SmartDataPtr<CalibData::CalibTest1>(m_pCalibDataSvc, CalibData::Test_Gen);
  
  std::string fullDeadPath = "/Calib/TKR_DeadChan/vanilla";
  DataObject *pDeadObject;
  

  m_pCalibDataSvc->retrieveObject(fullDeadPath, pDeadObject);

  CalibData::BadStrips* pDead = 0;
  pDead = dynamic_cast<CalibData::BadStrips *> (pDeadObject);
  if (!pDead) {
    log << MSG::ERROR << "Dynamic cast to BadStrips failed" << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNo = pDead->getSerNo();
  if (newSerNo != m_serDead) {
    log << MSG::INFO << "Processing new dead strips after retrieveObject" 
        << endreq;
    m_serDead = newSerNo;
    processNew(pDead, fullDeadPath);
  }
  m_pCalibDataSvc->updateObject(pDeadObject);

  pDead = 0;
  try {
    pDead = dynamic_cast<CalibData::BadStrips *> (pDeadObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to BadStrips after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pDead->getSerNo();
  if (newSerNo != m_serDead) {
    log << MSG::INFO << "Processing new dead strips after update" 
        << endreq;
    m_serDead = newSerNo;
    processNew(pDead, fullDeadPath);
  }


  // Now same for hot
  std::string fullHotPath = "/Calib/TKR_HotChan/vanilla";
  DataObject *pHotObject;

  m_pCalibDataSvc->retrieveObject(fullHotPath, pHotObject);

  CalibData::BadStrips* pHot = 0;
  pHot = dynamic_cast<CalibData::BadStrips *> (pHotObject);
  if (!pHot) {
    log << MSG::ERROR << "Dynamic cast to BadStrips failed" << endreq;
    return StatusCode::FAILURE;
  }

  newSerNo = pHot->getSerNo();
  if (newSerNo != m_serHot) {
    log << MSG::INFO << "Processing new hot strips after retrieveObject" 
        << endreq;
    m_serHot = newSerNo;
    processNew(pHot, fullHotPath);
  }
  m_pCalibDataSvc->updateObject(pHotObject);

  pHot = 0;
  try {
    pHot = dynamic_cast<CalibData::BadStrips *> (pHotObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to BadStrips after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pHot->getSerNo();
  if (newSerNo != m_serHot) {
    log << MSG::INFO << "Processing new hot strips after update" 
        << endreq;
    m_serHot = newSerNo;
    processNew(pHot, fullHotPath);
  }

  return StatusCode::SUCCESS;
}

void UseBadStrips::processNew(CalibData::BadStrips* pNew, 
                              const std::string& path) {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Retrieved with path " << path 
      << "Serial #" <<  pNew->getSerNo() << endreq; 
  log << MSG::INFO << "Vstart: " <<  (pNew->validSince()).hours()
      << "  Vend: " << (pNew->validTill()).hours() << endreq;
  log << MSG::INFO << "Bad type: " << pNew->getBadType() 
      << " has " << pNew->getBadTowerCount() << " bad towers " << endreq;

  pNew->traverse(m_visitor);
  
}

StatusCode UseBadStrips::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "          Finalize UseBadStrips "
      << endreq;
  
  return StatusCode::SUCCESS;
}

CalibData::eVisitorRet BadVisitor::badTower(unsigned int row, unsigned int col,
                                            int badness) {
  (*m_log) << MSG::INFO << "BadVisitor::badTower called with args "
           << "row = " << row << "col = "<< col 
           << "badness = " << badness << endreq;
  return CalibData::CONT;
}

CalibData::eVisitorRet BadVisitor::badPlane(unsigned int row, 
                                            unsigned int col, 
                                            unsigned int tray, bool top,
                                            int badness, 
                                            const CalibData::StripCol& strips)
{
  (*m_log) << MSG::INFO << "BadVisitor::badPlane called with "
           << "row = " << row << ", col = " << col << ", tray = "
           << tray << endreq;

  (*m_log) << MSG::INFO << "top = " << top << ", badness = " 
           << badness << endreq;
  (*m_log) << MSG::INFO << "Strip collection contains " << strips.size()
           << " strips. " << endreq;
  return CalibData::CONT;
}

//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibData/CalibModel.h"
#include "CalibData/Cal/CalCalibIntNonlin.h"
#include "CalibData/DacCol.h"
#include "idents/CalXtalId.h"   

/**
   @file UseIntNonlin.cxx
   Simple algorithm to test functioning of "the other" TDS, Cal integral
   nonlinearity data
*/



  /** 
   @class UseIntNonlin

   Algorithm exemplifying retrieval and use of Calorimeter int. nonlin. calib.
*/
class UseIntNonlin : public Algorithm {


public:
  UseIntNonlin(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper function called by execute
  void processNew(CalibData::CalCalibIntNonlin* pNew, const std::string& path);

  IDataProviderSvc* m_pCalibDataSvc;
  int               m_ser;
};


/// Instantiation of a static factory to create instances of this algorithm
static const AlgFactory<UseIntNonlin> Factory;
const IAlgFactory& UseIntNonlinFactory = Factory;


UseIntNonlin::UseIntNonlin(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
  m_ser(-1)
{
  // Declare properties here.

}


StatusCode UseIntNonlin::initialize() {
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


StatusCode UseIntNonlin::execute( ) {

  MsgStream log(msgSvc(), name());

  //  SmartDataPtr<CalibData::CalibTest1> test1(m_pCalibDataSvc,
  //                                        CalibData::Test_Gen);
  //  CalibData::CalibTest1* test1 = 
  //    SmartDataPtr<CalibData::CalibTest1>(m_pCalibDataSvc, CalibData::Test_Gen);
  
  std::string fullPath = "/Calib/CAL_IntNonlin/test";
  DataObject *pObject;
  

  m_pCalibDataSvc->retrieveObject(fullPath, pObject);

  CalibData::CalCalibIntNonlin* pIntNonlin = 0;
  pIntNonlin = dynamic_cast<CalibData::CalCalibIntNonlin *> (pObject);
  if (!pIntNonlin) {
    log << MSG::ERROR << "Dynamic cast to CalCalibIntNonlin failed" << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNo = pIntNonlin->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO 
        << "Processing new integral nonlinearity calibration after retrieveObject" 
        << endreq;
    m_ser = newSerNo;
    processNew(pIntNonlin, fullPath);
  }
  m_pCalibDataSvc->updateObject(pObject);

  pIntNonlin = 0;
  try {
    pIntNonlin = dynamic_cast<CalibData::CalCalibIntNonlin *> (pObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to CalCalibIntNonlin after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pIntNonlin->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new int. nonlin. after update" 
        << endreq;
    m_ser = newSerNo;
    processNew(pIntNonlin, fullPath);
  }

  return StatusCode::SUCCESS;
}


void UseIntNonlin::processNew(CalibData::CalCalibIntNonlin* pNew, 
                              const std::string& path) {
  using idents::CalXtalId;
  using CalibData::IntNonlin;
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
    //    unsigned range = 2;
    unsigned range = idents::CalXtalId::HEX8;
    unsigned face = 0;
    CalXtalId id(iTower, iLayer, iXtal);
    
    CalibData::RangeBase* pRange = pNew->getRange(id, range, face);
    
    IntNonlin* pIntNonlin = dynamic_cast<IntNonlin * >(pRange);
    log << MSG::INFO << "For tower = " << iTower << " layer = " << iLayer
        << " xtal = " << iXtal << endreq;
    log << MSG::INFO << "    range = " << range 
        << " face = " << face << endreq;

    const std::vector<float>* vals = pIntNonlin->getValues();
    unsigned nVal = vals->size();
    log << MSG::INFO << "Retrieved " << nVal
        << " integral nonlinearity values:"  << std::endl;
    for (unsigned iVal = 0; iVal < nVal; iVal++) {
      log << MSG::INFO << (*vals)[iVal] << " ";
    }
    log << MSG::INFO << std::endl;
    // log << MSG::INFO << "Averaged ped = " << pIntNonlin->getAvr() << endreq;
    log << MSG::INFO << "       error = " << pIntNonlin->getError() << endreq;

    /* Fetch dacs */

    for (unsigned iRange = 0; iRange < 4; iRange++) {
      CalibData::DacCol* dacCol = pNew->getDacCol(iRange);
      unsigned nD = (dacCol->getDacs())->size();
      log << MSG::INFO << "For range = " << iRange << " retrieved " 
          << nD << " dac values: " << endreq;

      for (unsigned iD = 0; iD < nD; iD++) {
        log << MSG::INFO << (*dacCol->getDacs())[iD] << " ";
      }
      log << endreq;
    }
    /*      Try another tower */
    iTower++;
    id = CalXtalId(iTower, iLayer, iXtal);
    
    pRange = pNew->getRange(id, range, face);

    if (pRange) {
      pIntNonlin = dynamic_cast<IntNonlin * >(pRange);
      log << MSG::INFO << "For tower = " << iTower << " layer = " << iLayer
          << " xtal = " << iXtal << endreq;
      log << MSG::INFO << "    range = " << range 
          << " face = " << face << endreq;
      
      // log << MSG::INFO << "Averaged ped = " 
      //     << pIntNonlin->getAvr() << endreq;
      log << MSG::INFO << "    error = " << pIntNonlin->getError() << endreq;
    }
  }
}

StatusCode UseIntNonlin::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "          Finalize UseIntNonlin "
      << endreq;
  
  return StatusCode::SUCCESS;
}


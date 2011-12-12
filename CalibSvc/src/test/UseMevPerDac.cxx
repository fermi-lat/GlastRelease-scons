//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibSvc/ICalibPathSvc.h"
#include "CalibData/Cal/CalMevPerDac.h"
#include "CalibData/Cal/Xpos.h"
#include "idents/CalXtalId.h"   

/**
   @file UseMevPerDac.cxx
   Simple algorithm to test functioning of "the other" TDS, Cal new 
   (Sept. 04) mevPerDac calibration data.
*/

  /** 
   @class UseMevPerDac

   Algorithm exemplifying retrieval and use of CAL_MevPerDac calibration.
*/
class UseMevPerDac : public Algorithm {

public:
  UseMevPerDac(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper function called by execute
  void processNew(CalibData::CalMevPerDacCol* pNew, const std::string& path);

  void printValSigs(MsgStream* log, const std::vector<CalibData::ValSig>* v, 
                    std::string title);

  void printValSig(MsgStream* log, const CalibData::ValSig* v, 
                   std::string title);
  IDataProviderSvc* m_pCalibDataSvc;
  ICalibPathSvc*    m_pCalibPathSvc;
  std::string       m_path;
  int               m_ser;
};


/// Instantiation of a static factory to create instances of this algorithm
//static const AlgFactory<UseMevPerDac> Factory;
//const IAlgFactory& UseMevPerDacFactory = Factory;
DECLARE_ALGORITHM_FACTORY(UseMevPerDac);

UseMevPerDac::UseMevPerDac(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
  m_ser(-1)
{
  // Declare properties here.

}


StatusCode UseMevPerDac::initialize() {
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
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_CAL_MevPerDac,
                                  std::string("test") );

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;
}


StatusCode UseMevPerDac::execute( ) {

  MsgStream log(msgSvc(), name());

  //  std::string fullPath = "/Calib/CAL_MevPerDac/test";
  DataObject *pObject;
  

  m_pCalibDataSvc->retrieveObject(m_path, pObject);

  CalibData::CalMevPerDacCol* pMevPerDac = 0;
  pMevPerDac = dynamic_cast<CalibData::CalMevPerDacCol *> (pObject);
  if (!pMevPerDac) {
    log << MSG::ERROR << "Dynamic cast to CalMevPerDacCol failed" << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNo = pMevPerDac->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO 
        << "Processing new MevPerDac calibration after retrieveObject" 
        << endreq;
    m_ser = newSerNo;
    processNew(pMevPerDac, m_path);
  }
  m_pCalibDataSvc->updateObject(pObject);

  pMevPerDac = 0;
  try {
    pMevPerDac = dynamic_cast<CalibData::CalMevPerDacCol *> (pObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to CalMevPerDacCol after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pMevPerDac->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new MevPerDac after update" 
        << endreq;
    m_ser = newSerNo;
    processNew(pMevPerDac, m_path);
  }

  return StatusCode::SUCCESS;
}


void UseMevPerDac::processNew(CalibData::CalMevPerDacCol* pNew, 
                              const std::string& path) {
  using idents::CalXtalId;
  using CalibData::CalMevPerDac;
  using CalibData::ValSig;
  using CalibData::Xpos;

  bool  done = false;

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Retrieved with path " << path << endreq
      << "Serial #" <<  pNew->getSerNo() << endreq; 
  log << MSG::INFO << "Vstart: " <<  (pNew->validSince()).hour(true)
      << "  Vend: " << (pNew->validTill()).hour(true) << endreq;
  
  if (!done) {
    done = true;
    short iTower = 0;
    short iLayer = 0;
    short iXtal = 2;
    //   range not applicable to this calibration type
    unsigned range = 0;
    unsigned face = 0;
    CalXtalId id(iTower, iLayer, iXtal);
    
    CalibData::RangeBase* pRange = pNew->getRange(id, range, face);
    
    CalMevPerDac* pMevPerDac = dynamic_cast<CalMevPerDac * >(pRange);
    log << MSG::INFO << "For tower = " << iTower << " layer = " << iLayer
        << " xtal = " << iXtal << endreq;

    if (!pMevPerDac) {
      log << MSG::INFO << "no calibration data found for this channel" 
          << endreq;
      return;
    }
    const std::vector<float>* positions = pNew->getXpos()->getVals();
    unsigned nVal = positions->size();

    unsigned nD = positions->size();
    log << MSG::INFO << "Xpos values are: " << endreq;    
    for (unsigned iD = 0; iD < nD; iD++) {
      log << MSG::INFO << (*positions)[iD] << " ";
    }
    log << endreq;

    log << MSG::INFO << "Retrieving " << nVal
        << " MevPerDac values and uncertainties:"  << endreq;

    const ValSig* pBig = pMevPerDac->getBig();
    const ValSig* pSm = pMevPerDac->getSmall();

    printValSig(&log, pBig, "Big: ");
    printValSig(&log, pSm, "Small: ");
    const std::vector<ValSig>* pRatio = 
      pMevPerDac->getBigSmallRatio(idents::CalXtalId::NEG);

    printValSigs(&log, pRatio, "bigSmallRatioN:");

    pRatio = 
      pMevPerDac->getBigSmallRatio(idents::CalXtalId::POS);

    printValSigs(&log, pRatio, "bigSmallRatioP:");

  }
}

void UseMevPerDac::printValSigs(MsgStream* log, 
                           const std::vector<CalibData::ValSig>* v, 
                           std::string title) {
  (*log) << MSG::INFO << "ValSig vector " << title << endreq;
  unsigned n = v->size();
  for (unsigned i = 0; i < n; i++) {
    (*log) << MSG::INFO << "val = " << (*v)[i].m_val 
           << "  sig = " << (*v)[i].m_sig << endreq;
  }
  return;
}

void UseMevPerDac::printValSig(MsgStream* log, 
                               const CalibData::ValSig* v, 
                               std::string title) {
  (*log) << MSG::INFO << "ValSig  " << title << endreq;
  (*log) << MSG::INFO << "val = " << v->m_val 
         << "  sig = " << v->m_sig << endreq;
  return;
}

StatusCode UseMevPerDac::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "Finalize UseMevPerDac "
      << endreq;
  
  return StatusCode::SUCCESS;
}


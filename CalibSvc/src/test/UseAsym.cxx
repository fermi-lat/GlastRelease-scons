//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibSvc/ICalibPathSvc.h"
#include "CalibData/Cal/CalAsym.h"
#include "CalibData/Cal/Xpos.h"
#include "idents/CalXtalId.h"   

/**
   @file UseAsym.cxx
   Simple algorithm to test functioning of "the other" TDS, Cal new 
   (Sept. 04) asymmetry calibration data.
*/

  /** 
   @class UseAsym

   Algorithm exemplifying retrieval and use of Calorimeter int. nonlin. calib.
*/
class UseAsym : public Algorithm {

public:
  UseAsym(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper function called by execute
  void processNew(CalibData::CalAsymCol* pNew, const std::string& path);

  void printValSigs(MsgStream* log, const std::vector<CalibData::ValSig>* v, 
                    std::string title);

  IDataProviderSvc* m_pCalibDataSvc;
  ICalibPathSvc*    m_pCalibPathSvc;
  std::string       m_path;
  int               m_ser;
};


/// Instantiation of a static factory to create instances of this algorithm
//static const AlgFactory<UseAsym> Factory;
//const IAlgFactory& UseAsymFactory = Factory;
DECLARE_ALGORITHM_FACTORY(UseAsym);


UseAsym::UseAsym(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
  m_ser(-1)
{
  // Declare properties here.

}


StatusCode UseAsym::initialize() {
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
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_CAL_Asym,
                                  std::string("test") );

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;
}


StatusCode UseAsym::execute( ) {

  MsgStream log(msgSvc(), name());

  //  std::string fullPath = "/Calib/CAL_Asym/test";
  DataObject *pObject;
  

  m_pCalibDataSvc->retrieveObject(m_path, pObject);

  CalibData::CalAsymCol* pAsym = 0;
  pAsym = dynamic_cast<CalibData::CalAsymCol *> (pObject);
  if (!pAsym) {
    log << MSG::ERROR << "Dynamic cast to CalAsymCol failed" << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNo = pAsym->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO 
        << "Processing new Asym calibration after retrieveObject" 
        << endreq;
    m_ser = newSerNo;
    processNew(pAsym, m_path);
  }
  m_pCalibDataSvc->updateObject(pObject);

  pAsym = 0;
  try {
    pAsym = dynamic_cast<CalibData::CalAsymCol *> (pObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to CalAsymCol after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pAsym->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new Asym after update" 
        << endreq;
    m_ser = newSerNo;
    processNew(pAsym, m_path);
  }

  return StatusCode::SUCCESS;
}


void UseAsym::processNew(CalibData::CalAsymCol* pNew, 
                              const std::string& path) {
  using idents::CalXtalId;
  using CalibData::CalAsym;
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

    // Range isn't applicable to Asym
    unsigned range = 0;
    unsigned face = 0;
    CalXtalId id(iTower, iLayer, iXtal);
    
    const std::vector<float>* positions = pNew->getXpos()->getVals();
    unsigned nVal = positions->size();
    
    unsigned nD = positions->size();
    log << MSG::INFO << "Xpos values are: " << endreq;    
    for (unsigned iD = 0; iD < nD; iD++) {
      log << MSG::INFO << (*positions)[iD] << " ";
    }
    log << endreq;
    
    log << MSG::INFO << "Retrieving " << nVal
        << " asym values and uncertainties:"  << std::endl;
    
    CalibData::RangeBase* pRange = pNew->getRange(id, range, face);
    
    CalAsym* pAsym = dynamic_cast<CalAsym * >(pRange);
    log << MSG::INFO << "For tower = " << iTower << " layer = " << iLayer
        << " xtal = " << iXtal << endreq;

    if (pAsym == 0) {
      log << MSG::INFO << "No calibration data for this channel" << endreq;
    }
    else {
    
      const std::vector<ValSig>* pV = pAsym->getBig();
      printValSigs(&log, pV, "Big:");
      
      pV = pAsym->getSmall();
      printValSigs(&log, pV, "Small:");
      
      pV = pAsym->getNSmallPBig();
      printValSigs(&log, pV, "NSmallPBig:");
      
      pV = pAsym->getPSmallNBig();
      printValSigs(&log, pV, "PSmallNBig:");
    }
  }
}

void UseAsym::printValSigs(MsgStream* log, 
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

StatusCode UseAsym::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "          Finalize UseAsym "
      << endreq;
  
  return StatusCode::SUCCESS;
}


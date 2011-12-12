//$Header$
#include <stdio.h>
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "CalibSvc/ICalibPathSvc.h"
#include "CalibData/Cal/CalCalibLightAsym.h"
#include "idents/CalXtalId.h"   

/**
   @file UseLightAsym.cxx
   Simple algorithm to test functioning of "the other" TDS, Cal light asym data

   NOTE:  this is currently (Sept. 2007) failing because there are no
   calibrations of type CAL_LightAsym in the database
*/

  /** 
   @class UseLightAsym

   Algorithm exemplifying retrieval and use of Calorimeter light asym calib.
*/
class UseLightAsym : public Algorithm {


public:
  UseLightAsym(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /// Helper function called by execute
  void processNew(CalibData::CalCalibLightAsym* pNew, const std::string& path);

  IDataProviderSvc* m_pCalibDataSvc;
  ICalibPathSvc*    m_pCalibPathSvc;
  std::string       m_path;
  int               m_ser;
};


/// Instantiation of a static factory to create instances of this algorithm
//static const AlgFactory<UseLightAsym> Factory;
//const IAlgFactory& UseLightAsymFactory = Factory;
DECLARE_ALGORITHM_FACTORY(UseLightAsym);


UseLightAsym::UseLightAsym(const std::string&  name, 
                 ISvcLocator*        pSvcLocator )
  : Algorithm(name, pSvcLocator), m_pCalibDataSvc(0), 
  m_ser(-1)
{
  // Declare properties here.

}


StatusCode UseLightAsym::initialize() {
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
    m_pCalibPathSvc->getCalibPath(ICalibPathSvc::Calib_CAL_LightAsym,
                                  std::string("vanilla") );

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  return StatusCode::SUCCESS;

}


StatusCode UseLightAsym::execute( ) {

  MsgStream log(msgSvc(), name());

  //  std::string fullPath = "/Calib/CAL_LightAsym/test";
  DataObject *pObject;
  

  m_pCalibDataSvc->retrieveObject(m_path, pObject);

  CalibData::CalCalibLightAsym* pLightAsym = 0;
  pLightAsym = dynamic_cast<CalibData::CalCalibLightAsym *> (pObject);
  if (!pLightAsym) {
    log << MSG::ERROR << "Dynamic cast to CalCalibLightAsym failed" << endreq;
    return StatusCode::FAILURE;
  }

  int newSerNo = pLightAsym->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO 
        << "Processing new light asym. calibration after retrieveObject" 
        << endreq;
    m_ser = newSerNo;
    processNew(pLightAsym, m_path);
  }
  m_pCalibDataSvc->updateObject(pObject);

  pLightAsym = 0;
  try {
    pLightAsym = dynamic_cast<CalibData::CalCalibLightAsym *> (pObject);
  }
  catch (...) {
    log << MSG::ERROR 
        << "Dynamic cast to CalCalibLightAsym after update failed" << endreq;
    return StatusCode::FAILURE;
  }
  newSerNo = pLightAsym->getSerNo();
  if (newSerNo != m_ser) {
    log << MSG::INFO << "Processing new pedestals after update" 
        << endreq;
    m_ser = newSerNo;
    processNew(pLightAsym, m_path);
  }

  return StatusCode::SUCCESS;
}


void UseLightAsym::processNew(CalibData::CalCalibLightAsym* pNew, 
                              const std::string& path) {
  using idents::CalXtalId;
  using CalibData::LightAsym;
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
    //    unsigned range = 2;
    unsigned range = idents::CalXtalId::HEX8;
    unsigned face = 0;
    CalXtalId id(iTower, iLayer, iXtal);
    
    CalibData::RangeBase* pRange = pNew->getRange(id, range, face);
    
    LightAsym* pLightAsym = dynamic_cast<LightAsym * >(pRange);
    log << MSG::INFO << "For tower = " << iTower << " layer = " << iLayer
        << " xtal = " << iXtal << endreq;
    log << MSG::INFO << "    range = " << range 
        << " face = " << face << endreq;

    const std::vector<float>* vals = pLightAsym->getValues();
    unsigned nVal = vals->size();
    log << MSG::INFO << "Retrieved " << nVal
        << " light asym values:"  << std::endl;
    for (unsigned iVal = 0; iVal < nVal; iVal++) {
      log << MSG::INFO << (*vals)[iVal] << " ";
    }
    log << MSG::INFO << std::endl;
    // log << MSG::INFO << "Averaged ped = " << pLightAsym->getAvr() << endreq;
    log << MSG::INFO << "       error = " << pLightAsym->getError() << endreq;

    /*      Try another tower */
    iTower++;
    id = CalXtalId(iTower, iLayer, iXtal);
    
    pRange = pNew->getRange(id, range, face);

    if (pRange) {
      pLightAsym = dynamic_cast<LightAsym * >(pRange);
      log << MSG::INFO << "For tower = " << iTower << " layer = " << iLayer
          << " xtal = " << iXtal << endreq;
      log << MSG::INFO << "    range = " << range 
          << " face = " << face << endreq;
      
      // log << MSG::INFO << "Averaged ped = " 
      //     << pLightAsym->getAvr() << endreq;
      log << MSG::INFO << "    error = " << pLightAsym->getError() << endreq;
    }
  }
}

StatusCode UseLightAsym::finalize( ) {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO 
      << "          Finalize UseLightAsym "
      << endreq;
  
  return StatusCode::SUCCESS;
}


// Implementation file for CalCalibSvc
//
// $Header $
//
// Author: Zach Fewtrell

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Incident.h"

#include "CalibData/DacCol.h"
#include "CalibData/CalibModel.h"
#include "CalCalibSvc.h"


static SvcFactory< CalCalibSvc > a_factory;
const ISvcFactory& CalCalibSvcFactory = a_factory; 

CalCalibSvc::CalCalibSvc(const std::string& name,ISvcLocator* Svc) 
  : Service(name,Svc) {


  m_intNonlinSerNo   = 0;

  // declare the properties
  declareProperty("CalibDataSvc",    m_calibDataSvcName = "CalibDataSvc");
  declareProperty("DefaultFlavor",   m_defaultFlavor    = "vanilla");
  declareProperty("FlavorGain",      m_flavorGain       = "");
  declareProperty("FlavorIntNonlin", m_flavorIntNonlin  = "");
  declareProperty("FlavorLightAsym", m_flavorLightAsym  = "");
  declareProperty("FlavorLightAtt",  m_flavorLightAtt   = "");
  declareProperty("FlavorMuSlope",   m_flavorMuSlope    = "");
  declareProperty("FlavorPed",       m_flavorPed        = "");
}

StatusCode  CalCalibSvc::queryInterface (const IID& riid, void **ppvIF) {
  if (IID_ICalCalibSvc == riid) {
	 *ppvIF = dynamic_cast<ICalCalibSvc*> (this);
	 return StatusCode::SUCCESS;
  }
  else {
	 return Service::queryInterface (riid, ppvIF);
  }
}

const IID&  CalCalibSvc::type () const {
  return IID_ICalCalibSvc;
}

StatusCode CalCalibSvc::initialize () 
{
  // Call super-class
  Service::initialize ();

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "initialize" << endreq;

  StatusCode  status = StatusCode::SUCCESS;

  
  // Grab pointer to CalibDataSvc
  status = service(m_calibDataSvcName, m_pCalibDataSvc, true);
  if ( !status.isSuccess() ) {
    log << MSG::ERROR << "Could not get CalibDataSvc" << endreq;
    return status;
  }

  // Grab interface pointers for CalibDataSvc's various interfaces
  // Query the IDataProvider interface of the calib data service
  status = m_pCalibDataSvc->queryInterface(IID_IDataProviderSvc, 
                                           (void**) &m_pDataProviderSvc);
  if ( !status.isSuccess() ) {
    log << MSG::ERROR 
        << "Could not query IDataProviderSvc interface of CalibDataSvc" 
        << endreq;
    return status;
  } else {
    log << MSG::DEBUG 
        << "Retrieved IDataProviderSvc interface of CalibDataSvc" 
        << endreq;
  }

  // Bind all of the properties for this service
  if ( (status = setProperties()).isFailure() ) {
	 log << MSG::ERROR << "Failed to set properties" << endreq;
  }

  // Post-init processing of jobOptions.
  // set default flavors unless otherwise specified
  if (!m_flavorGain.value().length())      m_flavorGain      = m_defaultFlavor;
  if (!m_flavorIntNonlin.value().length()) m_flavorIntNonlin = m_defaultFlavor;
  if (!m_flavorLightAsym.value().length()) m_flavorLightAsym = m_defaultFlavor;
  if (!m_flavorLightAtt.value().length())  m_flavorLightAtt  = m_defaultFlavor;
  if (!m_flavorMuSlope.value().length())   m_flavorMuSlope   = m_defaultFlavor;
  if (!m_flavorPed.value().length())       m_flavorPed       = m_defaultFlavor;
  
  log << MSG::DEBUG << "Initializing..."    << endreq;
  log << MSG::DEBUG << "  CalibDavaSvc   : " << m_calibDataSvcName << endreq;
  log << MSG::DEBUG << "  DefaultFlavor  : " << m_defaultFlavor    << endreq;
  log << MSG::DEBUG << "  FlavorGain     : " << m_flavorGain       << endreq;
  log << MSG::DEBUG << "  FlavorIntNonlin: " << m_flavorIntNonlin  << endreq;
  log << MSG::DEBUG << "  FlavorLightAsym: " << m_flavorLightAsym  << endreq;
  log << MSG::DEBUG << "  FlavorLightAtt : " << m_flavorLightAtt   << endreq;
  log << MSG::DEBUG << "  FlavorMuSlope  : " << m_flavorMuSlope    << endreq;
  log << MSG::DEBUG << "  FlavorPed      : " << m_flavorPed        << endreq;

  // Construct TDS data paths
  m_elecGainPath  = CalibData::CAL_ElecGain  + "/" + m_flavorGain.value();       
  m_intNonlinPath = CalibData::CAL_IntNonlin + "/" + m_flavorIntNonlin.value();  
  m_lightAsymPath = CalibData::CAL_LightAsym + "/" + m_flavorLightAsym.value();  
  m_lightAttPath  = CalibData::CAL_LightAtt  + "/" + m_flavorLightAtt.value();   
  m_muSlopePath   = CalibData::CAL_MuSlope   + "/" + m_flavorMuSlope.value();    
  m_pedPath       = CalibData::CAL_Ped       + "/" + m_flavorPed.value();        

  // Get ready to listen for BeginEvent
  IIncidentSvc* incSvc;
  status = service("IncidentSvc", incSvc, true);
  if (status.isSuccess() ) {
    int priority = 50;  // this should be lower priority (higher #?) than CalibDataSvc
    incSvc->addListener(this, "BeginEvent", priority);
  }
  else {
    log << MSG::ERROR << "Unable to find IncidentSvc" << endreq;
    return status;
  }

  return StatusCode::SUCCESS;
}

/// Inform that a new incident has occured
void CalCalibSvc::handle ( const Incident& inc ) { 
  MsgStream log(msgSvc(), name());

  if ((inc.type() == "BeginEvent")) {
    log << MSG::DEBUG << "New incident received" << endreq;
    log << MSG::DEBUG << "Incident source: " << inc.source() << endreq;
    log << MSG::DEBUG << "Incident type: " << inc.type() << endreq;
    
    // quick check if cache is valid & update if necessary.
    if (checkIntNonlinCache() != StatusCode::SUCCESS) {
      clearIntNonlinCache();
      initIntNonlinCache();
    }

  }
  return; 
}


StatusCode CalCalibSvc::finalize () {return StatusCode::SUCCESS;}

////////////////////////////////////////////////////////////////////////////////
/// BEGIN 'REAL' CODE - ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

StatusCode CalCalibSvc::getGain(const idents::CalXtalId &xtalId, 
										  idents::CalXtalId::XtalFace face,
										  idents::CalXtalId::AdcRange range, 
										  float &gain,
										  float &sig) {

  DataObject *pObject;
  CalibData::CalCalibGain *pGains  = 0;

  MsgStream log(msgSvc(), name());

  // Retrieve pointer to Gain tree from TDS
  if(m_pDataProviderSvc->retrieveObject(m_elecGainPath, pObject) == StatusCode::SUCCESS) {
	 pGains = dynamic_cast<CalibData::CalCalibGain *> (pObject);
	 if (!pGains) {
		log << MSG::ERROR << "Dynamic cast to CalCalibGain failed" << endreq;
		return StatusCode::FAILURE;
	 }
  } else{
	 log << MSG::ERROR << "Unable to retrieve gain from calib database" << endreq;
	 return StatusCode::FAILURE;
  }
  
  // Retrieve generic pointer to range specific data
  CalibData::RangeBase *pRangeBase;
  pRangeBase = pGains->getRange(xtalId, range,face);
  if (!pRangeBase) {
	 log << MSG::ERROR << "Unable to retrieve RangeBase for " << xtalId << endreq;
	 return StatusCode::FAILURE;
  }

  // recast for specific calibration type
  CalibData::Gain* pGain = dynamic_cast<CalibData::Gain *>(pRangeBase);
  if (!pGain) {
	 log << MSG::ERROR << "Dynamic cast to CalCalib::Gain failed" << endreq;
	 return StatusCode::FAILURE;
  }

  //values
  gain = pGain->getGain();
  sig  = pGain->getSig();

  return StatusCode::SUCCESS;
}

StatusCode CalCalibSvc::getIntNonlin(const idents::CalXtalId &xtalId, 
												 idents::CalXtalId::XtalFace face,
												 idents::CalXtalId::AdcRange range, 
												 const std::vector< float > *&vals,
												 const std::vector< unsigned > *&dacs,
												 float &error) {
  
  DataObject *pObject;
  CalibData::CalCalibIntNonlin *pIntNonlins;

  MsgStream log(msgSvc(), name());

  // Retrieve pointer to IntNonlin tree from TDS
  if(m_pDataProviderSvc->retrieveObject(m_intNonlinPath, pObject) == StatusCode::SUCCESS) {
	 pIntNonlins = dynamic_cast<CalibData::CalCalibIntNonlin *> (pObject);
	 if (!pIntNonlins) {
		log << MSG::ERROR << "Dynamic cast to CalCalibIntNonlin failed" << endreq;
		return StatusCode::FAILURE;
	 }
  } else{
	 log << MSG::ERROR << "Unable to retrieve IntNonlin from calib database" << endreq;
	 return StatusCode::FAILURE;
  }

  // Retrieve generic pointer to range specific data
  CalibData::RangeBase *pRangeBase;
  pRangeBase = pIntNonlins->getRange(xtalId, range,face);
  if (!pRangeBase) {
	 log << MSG::ERROR << "Unable to retrieve RangeBase for " << xtalId << endreq;
	 return StatusCode::FAILURE;
  }

  // recast for specific calibration type
  CalibData::IntNonlin* pIntNonlin = dynamic_cast<CalibData::IntNonlin *>(pRangeBase);
  if (!pIntNonlin) {
	 log << MSG::ERROR << "Dynamic cast to CalCalib::IntNonlin failed" << endreq;
	 return StatusCode::FAILURE;
  }

  //get vector of values
  const std::vector<float> *pIntNonlinVec = pIntNonlin->getValues();
  if (!pIntNonlin) {
	 log << MSG::ERROR << "Unable to get vector for IntNonlin vals" << endreq;
	 return StatusCode::FAILURE;
  }

  //get collection of associated DAC vals
  CalibData::DacCol *pIntNonlinDacCol = pIntNonlins->getDacCol(range);
  const std::vector<unsigned> *pIntNonlinDacVec;
  if (pIntNonlinDacCol) {
	 pIntNonlinDacVec = pIntNonlinDacCol->getDacs();
	 log << MSG::DEBUG << "pIntNonlinDacCol found /w " << pIntNonlinDacVec->size() << " elements\n" << endreq;
  } else {
	 log << MSG::ERROR << "No pIntNonlinDacCol found.\n" << endreq;
  }
  
  // check that Dac & ADC vectors are same size.
  if (pIntNonlinVec->size() != pIntNonlinDacVec->size()) {
	 log << MSG::ERROR << "pIntNonlinVec->size() != pIntNonlinDacVec->size() " 
        << "data=" << pIntNonlinVec->size() << " " 
        << "dac=" << pIntNonlinDacVec->size() << "xtalid=" << xtalId 
        << " face=" << face
        << " rng=" << range
        << endreq;
  }

  //whew, we made it this far... let's assign our outputs & leave!
  vals = pIntNonlinVec;
  dacs = pIntNonlinDacVec;
  error = pIntNonlin->getError();

  return StatusCode::SUCCESS;
}

StatusCode CalCalibSvc::getIntNonlin(const idents::CalXtalId &xtalId,
                                     idents::CalXtalId::XtalFace face,
                                     idents::CalXtalId::AdcRange range,
                                     const TSpline3 *&intNonlinSpline) {
  return retrieveIntNonlinSpline(xtalId,
                                 face,
                                 range,
                                 intNonlinSpline);
}
  
StatusCode CalCalibSvc::getLightAsym(const idents::CalXtalId &xtalId, 
												 idents::CalXtalId::XtalFace face,
												 idents::CalXtalId::AdcRange range, 
												 const std::vector< float > *&vals,
												 float &error) {
  DataObject *pObject;
  CalibData::CalCalibLightAsym *pLightAsyms  = 0;

  MsgStream log(msgSvc(), name());

  // Retrieve pointer to LightAsym tree from TDS
  if(m_pDataProviderSvc->retrieveObject(m_lightAsymPath, pObject) == StatusCode::SUCCESS) {
	 pLightAsyms = dynamic_cast<CalibData::CalCalibLightAsym *> (pObject);
	 if (!pLightAsyms) {
		log << MSG::ERROR << "Dynamic cast to CalCalibLightAsym failed" << endreq;
		return StatusCode::FAILURE;
	 }
  } else{
	 log << MSG::ERROR << "Unable to retrieve LightAsym from calib database" << endreq;
	 return StatusCode::FAILURE;
  }
  
  // Retrieve generic pointer to range specific data
  CalibData::RangeBase *pRangeBase;
  pRangeBase = pLightAsyms->getRange(xtalId, range,face);
  if (!pRangeBase) {
	 log << MSG::ERROR << "Unable to retrieve RangeBase for " << xtalId << endreq;
	 return StatusCode::FAILURE;
  }

  // recast for specific calibration type
  CalibData::LightAsym* pLightAsym = dynamic_cast<CalibData::LightAsym *>(pRangeBase);
  if (!pLightAsym) {
	 log << MSG::ERROR << "Dynamic cast to CalCalib::LightAsym failed" << endreq;
	 return StatusCode::FAILURE;
  }

  //get vector of values
  const std::vector<float> *pLightAsymVec = pLightAsym->getValues();
  if (!pLightAsym) {
	 log << MSG::ERROR << "Unable to get vector for LightAsym vals" << endreq;
	 return StatusCode::FAILURE;
  }

  //whew, we made it this far... let's assign our outputs & leave!
  vals = pLightAsymVec;
  error = pLightAsym->getError();

  return StatusCode::SUCCESS;
}
  
StatusCode CalCalibSvc::getLightAtt(const idents::CalXtalId &xtalId, 
												idents::CalXtalId::XtalFace face,
												idents::CalXtalId::AdcRange range, 
												float &att,
												float &norm) { 
  DataObject *pObject;
  CalibData::CalCalibLightAtt *pLightAtts  = 0;
  MsgStream log(msgSvc(), name());


  // Retrieve pointer to LightAtt tree from TDS
  if(m_pDataProviderSvc->retrieveObject(m_lightAttPath, pObject) == StatusCode::SUCCESS) {
	 pLightAtts = dynamic_cast<CalibData::CalCalibLightAtt *> (pObject);
	 if (!pLightAtts) {
		log << MSG::ERROR << "Dynamic cast to CalCalibLightAtt failed" << endreq;
		return StatusCode::FAILURE;
	 }
  } else{
	 log << MSG::ERROR << "Unable to retrieve LightAtt from calib database" << endreq;
	 return StatusCode::FAILURE;
  }
  
  // Retrieve generic pointer to range specific data
  CalibData::RangeBase *pRangeBase;
  pRangeBase = pLightAtts->getRange(xtalId, range,face);
  if (!pRangeBase) {
	 log << MSG::ERROR << "Unable to retrieve RangeBase for " << xtalId << endreq;
	 return StatusCode::FAILURE;
  }

  // recast for specific calibration type
  CalibData::LightAtt* pLightAtt = dynamic_cast<CalibData::LightAtt *>(pRangeBase);
  if (!pLightAtt) {
	 log << MSG::ERROR << "Dynamic cast to CalCalib::LightAtt failed" << endreq;
	 return StatusCode::FAILURE;
  }

  //whew, we made it this far... let's assign our outputs & leave!
  att = pLightAtt->getAtt();
  norm = pLightAtt->getNorm();

  return StatusCode::SUCCESS;
}

StatusCode CalCalibSvc::getMuSlope(const idents::CalXtalId &xtalId, 
											  idents::CalXtalId::XtalFace face,
											  idents::CalXtalId::AdcRange range, 
											  float &slope,
											  float &error) {
  DataObject *pObject;
  CalibData::CalCalibMuSlope *pMuSlopes  = 0;

  MsgStream log(msgSvc(), name());

  // Retrieve pointer to MuSlope tree from TDS
  if(m_pDataProviderSvc->retrieveObject(m_muSlopePath, pObject) == StatusCode::SUCCESS) {
	 pMuSlopes = dynamic_cast<CalibData::CalCalibMuSlope *> (pObject);
	 if (!pMuSlopes) {
		log << MSG::ERROR << "Dynamic cast to CalCalibMuSlope failed" << endreq;
		return StatusCode::FAILURE;
	 }
  } else{
	 log << MSG::ERROR << "Unable to retrieve MuSlope from calib database" << endreq;
	 return StatusCode::FAILURE;
  }
  
  // Retrieve generic pointer to range specific data
  CalibData::RangeBase *pRangeBase;
  pRangeBase = pMuSlopes->getRange(xtalId, range,face);
  if (!pRangeBase) {
	 log << MSG::ERROR << "Unable to retrieve RangeBase for " << xtalId << endreq;
	 return StatusCode::FAILURE;
  }

  // recast for specific calibration type
  CalibData::MuSlope* pMuSlope = dynamic_cast<CalibData::MuSlope *>(pRangeBase);
  if (!pMuSlope) {
	 log << MSG::ERROR << "Dynamic cast to CalCalib::MuSlope failed" << endreq;
	 return StatusCode::FAILURE;
  }

  //values
  slope = pMuSlope->getSlope();
  error  = pMuSlope->getError();

  return StatusCode::SUCCESS;
}
  
StatusCode CalCalibSvc::getPed(const idents::CalXtalId &xtalId, 
										 idents::CalXtalId::XtalFace face,
										 idents::CalXtalId::AdcRange range, 
										 float &avr,
										 float &sig,
										 float &cosAngle) {

  DataObject *pObject;
  CalibData::CalCalibPed *pPeds  = 0;

  MsgStream log(msgSvc(), name());

  // Retrieve pointer to Ped tree from TDS
  if(m_pDataProviderSvc->retrieveObject(m_pedPath, pObject) == StatusCode::SUCCESS) {
	 pPeds = dynamic_cast<CalibData::CalCalibPed *> (pObject);
	 if (!pPeds) {
		log << MSG::ERROR << "Dynamic cast to CalCalibPed failed" << endreq;
		return StatusCode::FAILURE;
	 }
  } else{
	 log << MSG::ERROR << "Unable to retrieve Peds from calib database" << endreq;
	 return StatusCode::FAILURE;
  }
  
  // Retrieve generic pointer to range specific data
  CalibData::RangeBase *pRangeBase;
  pRangeBase = pPeds->getRange(xtalId, range,face);
  if (!pRangeBase) {
	 log << MSG::ERROR << "Unable to retrieve RangeBase for " << xtalId << endreq;
	 return StatusCode::FAILURE;
  }

  // recast for specific calibration type
  CalibData::Ped* pPed = dynamic_cast<CalibData::Ped *>(pRangeBase);
  if (!pPed) {
	 log << MSG::ERROR << "Dynamic cast to CalCalib::Ped failed" << endreq;
	 return StatusCode::FAILURE;
  }

  //values
  avr      = pPed->getAvr();
  sig      = pPed->getSig();
  cosAngle = pPed->getCosAngle();
  return StatusCode::SUCCESS;
}

StatusCode CalCalibSvc::initIntNonlinCache() {
  DataObject *pObject;
  CalibData::CalCalibIntNonlin *pIntNonlins;

  MsgStream log(msgSvc(), name());

  // Retrieve pointer to IntNonlin tree from TDS
  if(m_pDataProviderSvc->retrieveObject(m_intNonlinPath, pObject) == StatusCode::SUCCESS) {
	 pIntNonlins = dynamic_cast<CalibData::CalCalibIntNonlin *> (pObject);
	 if (!pIntNonlins) {
		log << MSG::ERROR << "Dynamic cast to CalCalibIntNonlin failed" << endreq;
		return StatusCode::FAILURE;
	 }
  } else{
	 log << MSG::ERROR << "Unable to retrieve IntNonlin from calib database" << endreq;
	 return StatusCode::FAILURE;
  }

  // set current pointer
  m_intNonlinSerNo = pIntNonlins->getSerNo();
    
  return StatusCode::SUCCESS;
}

StatusCode CalCalibSvc::checkIntNonlinCache() {
  DataObject *pObject;
  CalibData::CalCalibIntNonlin *pIntNonlins;

  MsgStream log(msgSvc(), name());

  // Retrieve pointer to IntNonlin tree from TDS
  if(m_pDataProviderSvc->retrieveObject(m_intNonlinPath, pObject) == StatusCode::SUCCESS) {
	 pIntNonlins = dynamic_cast<CalibData::CalCalibIntNonlin *> (pObject);
	 if (!pIntNonlins) {
		log << MSG::ERROR << "Dynamic cast to CalCalibIntNonlin failed" << endreq;
		return StatusCode::FAILURE;
	 }
  } else{
	 log << MSG::ERROR << "Unable to retrieve IntNonlin from calib database" << endreq;
	 return StatusCode::FAILURE;
  }

  // check serial # for our cached data against current
  if (pIntNonlins->getSerNo() != m_intNonlinSerNo) {
    log << MSG::DEBUG << "IntNonlin data no longer valid" << endreq;
    return StatusCode::FAILURE;
  }
  
  return StatusCode::SUCCESS;
}

StatusCode CalCalibSvc::clearIntNonlinCache() {
  m_intNonlinList.Delete();
  m_intNonlinSerNo = 0;
  
  return StatusCode::SUCCESS;
}

StatusCode CalCalibSvc::loadIntNonlinSpline(const idents::CalXtalId &xtalId,
                                            idents::CalXtalId::XtalFace face,
                                            idents::CalXtalId::AdcRange range) {

  MsgStream log(msgSvc(), name());

  std::string grName;
  generateSplineName(grName,                               
                     std::string("CAL_IntNonlin"),
                     xtalId, face, range);
  
  // see if we already have a TSpline3 object and delete it.
  TSpline3 *pSpline = (TSpline3 *)m_intNonlinList.FindObject(grName.c_str());
  if (pSpline) pSpline->Delete();
  delete pSpline;
  pSpline = 0;


  // now create new TSpline3y object
  const std::vector<float> *pIntNonlinVec;
  const std::vector<unsigned> *pIntNonlinDacVec;
  float error;

  if (getIntNonlin(xtalId, face, range, 
                   pIntNonlinVec, pIntNonlinDacVec, error)
      != StatusCode::SUCCESS)
    return StatusCode::FAILURE;


  // create temporary array for storing doubles for input into tspline.
  int n = pIntNonlinVec->size();
  double *pIntNonlin = new double[n];
  double *pDac = new double[pIntNonlinDacVec->size()];
                                  
  for (int i = 0; i < n; i++) {
    pIntNonlin[i] = (*pIntNonlinVec)[i];
    pDac[i] = (*pIntNonlinDacVec)[i];
  }

  // create spline object.
  TSpline3 *mySpline = new TSpline3(grName.c_str(), pDac, pIntNonlin, n);
  delete pIntNonlin;
  delete pDac;

  mySpline->SetName(grName.c_str());
  m_intNonlinList.Add(mySpline);

  return StatusCode::SUCCESS;
}

StatusCode CalCalibSvc::retrieveIntNonlinSpline(const idents::CalXtalId &xtalId,
                                                idents::CalXtalId::XtalFace face,
                                                idents::CalXtalId::AdcRange range,
                                                const TSpline3 *&pSpline) {
  MsgStream log(msgSvc(), name());

  std::string grName; 
  generateSplineName(grName,                               
                     std::string("CAL_IntNonlin"),
                     xtalId, face, range);
  
  pSpline = (TSpline3 *)m_intNonlinList.FindObject(grName.c_str());
  if (pSpline != NULL) return StatusCode::SUCCESS;

  // object does not exist in cache
  if (loadIntNonlinSpline(xtalId,face,range) != StatusCode::SUCCESS) return StatusCode::FAILURE;
  pSpline = (TSpline3 *)m_intNonlinList.FindObject(grName.c_str());
  
  return StatusCode::SUCCESS;
}


std::string &CalCalibSvc::generateSplineName(std::string &grName,              // output
                                             std::string grType, 
                                             const idents::CalXtalId &xtalId,
                                             idents::CalXtalId::XtalFace face,
                                             idents::CalXtalId::AdcRange range) {
  std::ostringstream tmp_stream;
  
  tmp_stream << grType 
             << 'x' << xtalId
             << 'f' << face
             << 'r' << range
             << std::ends;

  grName = std::string(tmp_stream.str());

  return grName;
}


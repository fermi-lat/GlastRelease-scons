/*
@file TrgConfigSvc.cxx

@brief keeps track of the GEM trigger configuration
@author Martin Kocian

$Header$

*/

#include "./ConfigSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"

#include "GaudiKernel/SmartDataPtr.h"
#include "LdfEvent/LsfMetaEvent.h"
//#include "configData/db/LatcDBImplOld.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "mootCore/MootQuery.h"
#include "CalibData/Moot/MootData.h"
#include "MootSvc/IMootSvc.h"
#include "configData/fsw/FswEfcSampler.h"
#include "configData/gem/TrgConfigParser.h"
#include "facilities/commonUtilities.h"

#include <set>
#include <stdlib.h>
#include <assert.h>


// declare the service factories for the ConfigSvc
static SvcFactory<ConfigSvc> a_factory;
const ISvcFactory& ConfigSvcFactory = a_factory; 


ConfigSvc::ConfigSvc(const std::string& name,ISvcLocator* svc) 
  : Service(name,svc),m_mootSvc(0),m_trgConfig(0),m_noMOOTMask(0),m_mootKey(0),m_latcKey(0){
  // Purpose and Method: Constructor - Declares and sets default properties
  //                     
  // Inputs: service name and locator 
  //         
  // Outputs: None
  // Dependencies: None
  // Restrictions and Caveats:  None
  
  // declare the properties  
  declareProperty("GammaFilterXml",  m_gammaFilterXmlFile="");
  declareProperty("DgnFilterXml",  m_dgnFilterXmlFile="");
  declareProperty("MipFilterXml",  m_mipFilterXmlFile="");
  declareProperty("HipFilterXml",  m_hipFilterXmlFile="");
  declareProperty("GemXml", m_gemXmlFile="");
  declareProperty("RoiXml", m_roiXmlFile="");
}

StatusCode  ConfigSvc::queryInterface (const InterfaceID& riid, void **ppvIF){
  if (IID_IConfigSvc == riid) {
    *ppvIF = dynamic_cast<IConfigSvc*> (this);
    return StatusCode::SUCCESS;
  } else {
    return Service::queryInterface (riid, ppvIF);
  }
}

const InterfaceID&  ConfigSvc::type () const {
  return IID_IConfigSvc;
}

StatusCode ConfigSvc::initialize () {
  // Purpose and Method: Initialize the lists of dead units
  //                     
  // Inputs: None        
  // Outputs: None
  // Dependencies: None
  // Restrictions and Caveats:  None
  
  StatusCode  sc = StatusCode::SUCCESS;
  
  Service::initialize();
  
  // Open the message log
  MsgStream log( msgSvc(), name() );

  // MOOT service
  if ((sc = service("MootSvc",m_mootSvc,true)).isFailure()) { 
    log << MSG::ERROR << "Failed to find MOOT service" << endreq;
    return StatusCode::FAILURE;
  }
    
  // Bind all of the properties for this service
  if ( (sc = setProperties()).isFailure() ) {
    log << MSG::ERROR << "Failed to set properties" << endreq;
    return StatusCode::FAILURE;
  }

  if ( m_gemXmlFile.value() != "" || m_roiXmlFile.value() != "" ) {
    m_noMOOTMask |= ( 1 << GEM );
    m_noMOOTMask |= ( 1 << ROI );    
    if ( ! readTrgConfig( m_gemXmlFile.value().c_str(), m_roiXmlFile.value()) ) {
      log << MSG::ERROR << "Failed to read GEM configuration from " << m_gemXmlFile.value() << endreq;
      return StatusCode::FAILURE;
    } else {
      log << MSG::WARNING << "Using GEM configuration from file " << m_gemXmlFile.value() << " instead of MOOT." << endreq;
    }
  } else {
    log << MSG::WARNING << "Using GEM configuration from MootSvc" << endreq;
  }
  if ( m_gammaFilterXmlFile.value() != "" ) {
    m_noMOOTMask |= ( 1 << enums::Lsf::GAMMA );   
    if ( ! readEfcFromFile( MOOT::LPA_MODE_ALL, enums::Lsf::GAMMA, m_gammaFilterXmlFile.value() ) ) {
      log << MSG::ERROR << "Failed to read Gamma filter configuration from " << m_gammaFilterXmlFile.value() << endreq;
      return StatusCode::FAILURE;
    } else {
      log << MSG::WARNING << "Using GAMMA filter configuration from file " << m_gammaFilterXmlFile.value() << " instead of MOOT." << endreq;
    }
  } else {
    log << MSG::WARNING << "Using GAMMA filter configuration from MootSvc" << endreq;
  }
  if ( m_dgnFilterXmlFile.value() != "" ) {
     m_noMOOTMask |= ( 1 << enums::Lsf::DGN );  
    if ( ! readEfcFromFile( MOOT::LPA_MODE_ALL, enums::Lsf::DGN, m_dgnFilterXmlFile.value() ) ) {
      log << MSG::ERROR << "Failed to read DGN filter configuration from " << m_dgnFilterXmlFile.value() << endreq;
      return StatusCode::FAILURE;
    } else {
      log << MSG::WARNING << "Using DGN filter configuration from file " << m_dgnFilterXmlFile.value() << " instead of MOOT." << endreq;
    }
  } else {
    log << MSG::WARNING << "Using DGN filter configuration from MootSvc" << endreq;
  }
  if ( m_mipFilterXmlFile.value() != "" ) {
    m_noMOOTMask |= ( 1 << enums::Lsf::MIP );  
    if ( ! readEfcFromFile( MOOT::LPA_MODE_ALL, enums::Lsf::MIP, m_mipFilterXmlFile.value() ) ) {
      log << MSG::ERROR << "Failed to read MIP filter configuration from " << m_mipFilterXmlFile.value() << endreq;
      return StatusCode::FAILURE;
    } else {
      log << MSG::WARNING << "Using MIP filter configuration from file " << m_mipFilterXmlFile.value() << " instead of MOOT." << endreq;
    }
  } else {
    log << MSG::WARNING << "Using MIP filter configuration from MootSvc" << endreq;
  }
  if ( m_hipFilterXmlFile.value() != "" ) {
    m_noMOOTMask |= ( 1 << enums::Lsf::HIP );  
    if ( ! readEfcFromFile( MOOT::LPA_MODE_ALL, enums::Lsf::HIP, m_hipFilterXmlFile.value() ) ) {
      log << MSG::ERROR << "Failed to read HIP filter configuration from " << m_hipFilterXmlFile.value() << endreq;
      return StatusCode::FAILURE;
    } else {
      log << MSG::WARNING << "Using HIP filter configuration from file " << m_hipFilterXmlFile.value() << " instead of MOOT." << endreq;
    }
  } else {
    log << MSG::WARNING << "Using HIP filter configuration from MootSvc" << endreq;
  }
  return sc;
}


unsigned ConfigSvc::getMootKey() const {
  newMootKey();
  return m_mootKey;
}


const TrgConfig* ConfigSvc::getTrgConfig() const {
  
  std::string gemFileName;
  std::string roiFileName;
  unsigned checkMOOT = ( (1 << GEM) | (1 << ROI) );
  if ( (checkMOOT &  m_noMOOTMask) == checkMOOT ) {
    // no need to check moot, just return the TrgConfig
    if ( m_trgConfig == 0 ) {
      gemFileName = m_gemXmlFile.value();
      roiFileName = m_roiXmlFile.value();
      if ( ! readTrgConfig( gemFileName, roiFileName ) ) {
	return 0;
      } 
    }
    return m_trgConfig;
  }
  if ( ! newMootKey() &&  m_trgConfig != 0) {
    // same key, just return the TrgConfig
    return m_trgConfig;
  }
  unsigned latcKey(0);
  if ( m_noMOOTMask & ( 1 << GEM ) ) {
    gemFileName = m_gemXmlFile.value();
  } else {
    if ( m_mootSvc->noMoot() ) return 0;
    const CalibData::MootParm* gem = m_mootSvc->getGemParm(latcKey);    
    //getFullPath( gem->getSrc(), gemFileName );
    gemFileName = gem->getSrc();
  }
  if ( m_noMOOTMask & ( 1 << ROI ) ) {
    roiFileName = m_roiXmlFile.value();
  } else {
    if ( m_mootSvc->noMoot() ) return 0;
    const CalibData::MootParm* roi = m_mootSvc->getRoiParm(latcKey);
    //getFullPath( roi->getSrc(), roiFileName );
    roiFileName = roi->getSrc();
  }
  if ( ! readTrgConfig( gemFileName, roiFileName ) ) {
    return 0;
  }
  return m_trgConfig;
}


const FswEfcSampler* ConfigSvc::getFSWPrescalerInfo(enums::Lsf::Mode mode, unsigned handlerId ) const {

  MOOT::LpaMode mootMode = MOOT::LPA_MODE_count;
  switch (mode) {
  case enums::Lsf::Normal: 
    mootMode = MOOT::LPA_MODE_NORMAL; break;
  case enums::Lsf::TOO:    
    mootMode = MOOT::LPA_MODE_TOO; break;
  case enums::Lsf::GRB0:
    mootMode = MOOT::LPA_MODE_GRB0; break;
  case enums::Lsf::GRB1:
    mootMode = MOOT::LPA_MODE_GRB1; break;
  case enums::Lsf::GRB2:
    mootMode = MOOT::LPA_MODE_GRB2; break;
  case enums::Lsf::Solar:
    mootMode = MOOT::LPA_MODE_SOLAR; break;
  case enums::Lsf::Calibration:
    mootMode = MOOT::LPA_MODE_CALIBRATION; break;
  case enums::Lsf::Diagnostic:    
    mootMode = MOOT::LPA_MODE_DIAGNOSTIC; break;
  default:
    break;
  }
  return getFswSampler(mootMode,handlerId);
}


void ConfigSvc::getFullPath(const std::string& mootPath, std::string& fullPath) const {
  std::string archEnv("MOOT_ARCHIVE");
  std::string localStr = facilities::commonUtilities::getEnvironment( archEnv );
  fullPath = facilities::commonUtilities::joinPath( facilities::commonUtilities::getEnvironment( archEnv ),
						    mootPath );
}



/// check to see if we have a new MOOT key
bool ConfigSvc::newMootKey() const {
  if ( m_mootSvc->noMoot() ) return false;
  unsigned newMootKey = m_mootSvc->getMootConfigKey();
  unsigned newLatcKey = m_mootSvc->getHardwareKey();
  if ( m_mootKey == newMootKey && m_latcKey == newLatcKey ) return false;
  m_mootKey = newMootKey;  
  m_latcKey = newLatcKey;
  clearCache();
  return true;
}


/// clear our the read data strutures
void ConfigSvc::clearCache() const {
  delete m_trgConfig;
  m_trgConfig = 0;
  std::set<FswEfcSampler*> deleted;
  for ( std::map<unsigned,FswEfcSampler*>::iterator itr = m_fswEfcCache.begin(); itr !=  m_fswEfcCache.end(); itr++ ) {
    FswEfcSampler* sampler = itr->second;
    if ( deleted.find(sampler) == deleted.end() ) {
      delete itr->second;
      deleted.insert(sampler);
    }
  }
  m_fswEfcCache.clear();
}



/// read the TrgConfig
bool ConfigSvc::readTrgConfig( const std::string& gemFileName, const std::string& roiFileName ) const {
  TrgConfigParser parser;
  m_trgConfig = new TrgConfig;
  if ( gemFileName != "" ) {
    if ( parser.parse(m_trgConfig,gemFileName.c_str()) != 0 ) {
      delete m_trgConfig;
      m_trgConfig = 0;
      return false;
    }
  }
  if ( roiFileName != "" ) {
    if ( parser.parse(m_trgConfig,roiFileName.c_str()) != 0 ) {
      delete m_trgConfig;
      m_trgConfig = 0;
      return false;
    }
  }
  return true;
}



/// function to get the prescale factor
FswEfcSampler* ConfigSvc::getFswSampler(unsigned lpaMode, unsigned handlerId) const {

  // Open the message log
  MsgStream log( msgSvc(), name() );

  FswEfcSampler* sampler(0);

  // Check to see if we are reading that handler from a file instead of from MOOT
  if ( m_noMOOTMask & ( 1 << handlerId ) ) {
    return getFswSamplerFromCache(lpaMode,handlerId);
  }
  // Check the cache
  if ( ! newMootKey() ) {
    // it is an old key, return the cached value
    sampler = getFswSamplerFromCache(lpaMode,handlerId);
    // No previous cache, try and get it.
    if ( sampler != 0 ) return sampler;
  } 
  if ( m_mootSvc->noMoot() ) {
    log << MSG::ERROR << "Need to get FswHandler from MOOT, but IMootSvc::noMoot() is set." << endreq;
    return 0;
  }
  // new value, read the xml file into the cache
  std::string handlerName; 
  CalibData::MootFilterCfg* filterCfg = m_mootSvc->getActiveFilter(lpaMode,handlerId,handlerName);
  if ( filterCfg == 0 ) {
    log << MSG::ERROR << "No active filter for mode " << lpaMode << " and handlerId " << handlerId << endreq;
    return 0;
  }
  std::string fullPath;
  getFullPath( filterCfg->getSrcPath(), fullPath );

  sampler = readEfcFromFile(lpaMode,handlerId,fullPath);
  if ( sampler == 0 ) {
    log << MSG::ERROR << "Failed to read sampler data from file " << fullPath << endreq;
  } else {
    log << MSG::INFO << "Read sampler data from file " << fullPath << endreq;
  }
  return sampler;
}


/// function to get the cached prescale factors
FswEfcSampler* ConfigSvc::getFswSamplerFromCache(unsigned lpaMode, unsigned handlerId) const {
  unsigned slotKey = filterSlotKey(lpaMode,handlerId);
  std::map<unsigned,FswEfcSampler*>::iterator itr = m_fswEfcCache.find(slotKey);
  FswEfcSampler* retVal = itr != m_fswEfcCache.end() ? itr->second : 0;
  return retVal;
}


/// read an xml file with the prescaler info 
FswEfcSampler* ConfigSvc::readEfcFromFile(unsigned lpaMode, unsigned handlerId, const std::string& fileName ) const {    
  FswEfcSampler* sampler = FswEfcSampler::makeFromXmlFile(fileName.c_str());
  if ( sampler == 0 ) return 0;
  if ( lpaMode == MOOT::LPA_MODE_ALL ) {
    for ( unsigned i =  MOOT::LPA_MODE_NORMAL; i < MOOT::LPA_MODE_count; i++ ) {
      unsigned slotKey = filterSlotKey(i,handlerId);
      std::map<unsigned,FswEfcSampler*>::iterator itr = m_fswEfcCache.find(slotKey);
      if ( itr == m_fswEfcCache.end() ) {
	m_fswEfcCache[slotKey] = sampler;
      } else {
	delete itr->second;
	itr->second = sampler;
      } 
    }
  } else {
    unsigned slotKey = filterSlotKey(lpaMode,handlerId);
    std::map<unsigned,FswEfcSampler*>::iterator itr = m_fswEfcCache.find(slotKey);
    if ( itr == m_fswEfcCache.end() ) {
      m_fswEfcCache[slotKey] = sampler;
    } else {
      delete itr->second;
      itr->second = sampler;
    } 
  }

  MsgStream log( msgSvc(), name() );
  log << MSG::INFO << "Read Filter Configuration for Mode: " << lpaMode << "; Handler: " 
      << handlerId << "; from file: " << fileName << std::endl
      << *sampler << endreq;
  
  return sampler;
}



	    
// handle "incidents"

StatusCode ConfigSvc::finalize( ) {
  MsgStream log(msgSvc(), name());  
  clearCache();
  return StatusCode::SUCCESS;
}

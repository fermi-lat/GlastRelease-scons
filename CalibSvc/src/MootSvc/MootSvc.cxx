//$Header$
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

#include <utility>   // for std::pair
// <vector>, <string> included by MootSvc.h

#include "MootSvc.h"
#include "mootCore/MoodConnection.h"
#include "mootCore/MootQuery.h"
#include "mootCore/FileDescrip.h"
#include "facilities/Util.h"
#include "facilities/commonUtilities.h"
#include "CalibSvc/MootTds.h"
#include "../CalibDataSvc/CalibCLIDNode.h"

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/GenericAddress.h"
#include "GaudiKernel/IConverter.h"

#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IValidity.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"

/// Instantiation of a static factory to create instances of this service
static SvcFactory<MootSvc>          MootSvc_factory;
const ISvcFactory& MootSvcFactory = MootSvc_factory;

// For associating enum with parm class names
/*
struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return std::strcmp(s1, s2) == 0;
  }
};

#ifdef _WIN32
typedef hash_map<const char*, int, 
                 hash<const char*, eqstr> > HashMap;
#else
typedef hash_map<const char*, int> HashMap;
                 //                 hash<const char*, eqstr> > HashMap;
#endif
*/
MootSvc::MootSvc( const std::string& name, ISvcLocator* svc)
  : ConversionSvc (name, svc, MOOT_StorageType), m_q(0), m_c(0), m_log(0)
{
  declareProperty("MootArchive", m_archive = std::string("") );
  m_latcPath.clear();
}

MootSvc::~MootSvc(){ }

StatusCode MootSvc::initialize()
{
  // Initialize base class
  StatusCode sc = ConversionSvc::initialize();
  if ( !sc.isSuccess() ) return sc;

  m_log = new MsgStream(msgSvc(), "MootSvc");

  (*m_log) << MSG::INFO << "Specific initialization starting" << endreq;

  // Locate the Calib Data Service.  Since it inherits from DataSvc
  // it has to implement IDataProviderSvc
  IDataProviderSvc* pCDS = 0;
  sc = serviceLocator()->getService 
    ("CalibDataSvc",  IID_IDataProviderSvc, (IInterface*&)pCDS);
  if ( !sc.isSuccess() ) {
    (*m_log) << MSG::ERROR << "Could not locate CalibDataSvc" << endreq;
    return sc;
  }

  // Set the CalibDataSvc as data provider service
  sc = setDataProvider(pCDS);
  if ( !sc.isSuccess() ) {
    (*m_log) << MSG::ERROR << "Could not set data provider" << endreq;
    return sc;
  }

  // Query the DataProviderSvc and  DataManagerSvc interfaces of the 
  // calib data service.  Don't think we need detDataSvc interface
  //  sc = pCDS->queryInterface(IID_IDetDataSvc, (void**) &m_detDataSvc);
  sc = pCDS->queryInterface(IID_IDataProviderSvc, (void**) &m_provider);
  if ( !sc.isSuccess() ) {
    (*m_log) << MSG::ERROR 
	<< "Cannot query IDataProviderSvc interface of CalibDataSvc" << endreq;
    return sc;
  } else {

    (*m_log) << MSG::DEBUG 
             << "Retrieved IDataProviderSvc interface of CalibDataSvc" 
	<< endreq;
  }

  sc = pCDS->queryInterface(IID_IDataManagerSvc, (void**) &m_dataManager);
  if ( !sc.isSuccess() ) {
    (*m_log) << MSG::ERROR 
	<< "Cannot query IDataManager interface of CalibDataSvc" << endreq;
    return sc;
  } else {

    (*m_log) << MSG::DEBUG 
             << "Retrieved IDataManagerSvc interface of CalibDataSvc" 
             << endreq;
  }

  // Not sure this is necessary, but leave it in for now
  sc = pCDS->queryInterface(IID_IInstrumentName, (void**) &m_instrSvc);
  if ( !sc.isSuccess() ) {
    (*m_log) << MSG::ERROR 
	<< "Cannot query IInstrumentName interface of CalibDataSvc" << endreq;
    return sc;
  } else {
    (*m_log) << MSG::DEBUG 
        << "Retrieved IInstrumentNameSvc interface of CalibDataSvc" << endreq;
  }

  // Locate IConversionSvc interface of the DetectorPersistencySvc
  sc = serviceLocator()->service 
    ("DetectorPersistencySvc", m_detPersSvc, true);
  if ( !sc.isSuccess() ) {
    (*m_log) << MSG::ERROR 
	<< "Cannot locate IConversionSvc interface of DetectorPersistencySvc"
	<< endreq;
    return sc;
  } else {
    (*m_log) << MSG::DEBUG 
	<< "Retrieved IConversionSvc interface of DetectorPersistencySvc"
	<< endreq;
  }
  IAddressCreator* addrCreator;

  // Query the IAddressCreator interface of the detector persistency service
  sc = m_detPersSvc->queryInterface(IID_IAddressCreator, 
				    (void**) &addrCreator);
  if ( !sc.isSuccess() ) {
    (*m_log) << MSG::ERROR 
	<< "Cannot access IAddressCreator interface of DetectorPersistencySvc" 
	<< endreq;
    return sc;
  } else {
    (*m_log) << MSG::DEBUG 
	<< "Retrieved IAddressCreator interface of DetectorPersistencySvc" 
	<< endreq;
  }
  //  log << MSG::DEBUG 
  //    << "Set it as the address creator of MootSvc" << endreq;
  sc = setAddressCreator(addrCreator);
  if ( !sc.isSuccess() ) {
    (*m_log) << MSG::ERROR 	<< "Cannot set the address creator" << endreq;
    return sc;
  }

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  if ( !sc.isSuccess() ) {
    (*m_log) << MSG::ERROR << "Could not set jobOptions properties" << endreq;
    return sc;
  }
  (*m_log) << MSG::DEBUG << "Properties were read from jobOptions" << endreq;


  /*
    Might want to add something analogous here to force use of
    user-supplied LATC master key
  if (!m_useEventTime) {  // special diagnostic mode
  }

  */
  // Make a MOOT::MootQuery instance
  // Conceivably, could start up a different conversion service, depending 
  // on job options parameters, which would look very much like this one
  // except for having a different way to access metadata.

  m_q = getConnection(true);    // for now, always verbose
  if (!m_q) return StatusCode::FAILURE;

  (*m_log) << MSG::INFO << "Specific initialization completed" << endreq;
  return sc;
}


StatusCode MootSvc::finalize()
{
  (*m_log) << MSG::DEBUG << "Finalizing" << endreq;
  if (m_q) delete m_q;
  if (m_c) delete m_c;
  delete m_log;
  m_q = 0;  m_c = 0; m_log = 0;
    
  return ConversionSvc::finalize();
}

MOOT::MootQuery* MootSvc::getConnection(bool verbose) {
  if (m_q) return m_q;

  const std::string slacDefault("/afs/slac/glast/g/moot/archive-mood/");

  std::string archEnv("$(MOOT_ARCHIVE)");
  std::string archEnvName("MOOT_ARCHIVE");

  bool envSet = false;

  // To make a connection, need definition for env. var MOOT_ARCHIVE
  //   If path has been supplied in job options, do setenv for it
  //     Special value of "*" means use default 
  //     (/afs/slac/g/glast/moot/archive-mood)
  //   else if $MOOT_ARCHIVE already has def, use that
  //   else try default value above 

  if (m_archive.size() == 0 ) {
    // Check to see if MOOT_ARCHIVE has a value.  If not, set m_archive to 
    int nExpand = facilities::Util::expandEnvVar(&archEnv);
    if (nExpand > 0) envSet = true;
    else m_archive = slacDefault;
  }
  if (!envSet) {
    if (m_archive == std::string("*"))   { // use slac default
      m_archive = slacDefault;
    }
    facilities::commonUtilities::setEnvironment(archEnvName, m_archive, true);
  }
  // Maybe should have some logic here to get verbose connection
  // depending on debug level?
  m_c = new MOOT::MoodConnection(false, verbose);
  if (m_c) m_q = new MOOT::MootQuery(m_c);

  if (!m_q) {
    (*m_log) << MSG::ERROR << "Could not open connection to MOOT dbs" << endreq;
  }
  
  return m_q;
}

unsigned MootSvc::getLatcParmMaxCnt() {
  if (m_latcParmMap) return m_latcParmMap->size();
  else return 0;
}

StatusCode MootSvc::getPrecincts() {
  
  m_latcParmMap = new HashMap();

  m_prNames.resize(PR_count);
  // We think we know what they should be called.  Initialize, then confirm
  // that all of them are real precinct names by checking against db.

  m_prNames[PR_generic] = std::string("generic");
  m_prNames[PR_ACD_Mode] = std::string("ACD_Mode");
  m_prNames[PR_ACD_Bias] = std::string("ACD_Bias");
  m_prNames[PR_ACD_Hld] = std::string("ACD_Hld");
  m_prNames[PR_ACD_PHA] = std::string("ACD_PHA");
  m_prNames[PR_ACD_Veto] = std::string("ACD_Veto");
  m_prNames[PR_ACD_Timing] = std::string("ACD_Timing");
  m_prNames[PR_CAL_Mode] = std::string("CAL_Mode");
  m_prNames[PR_CAL_Timing] = std::string("CAL_Timing");
  m_prNames[PR_CAL_LAC] = std::string("CAL_LAC");
  m_prNames[PR_CAL_FLE] = std::string("CAL_FLE");
  m_prNames[PR_CAL_FHE] = std::string("CAL_FHE");
  m_prNames[PR_CAL_ULD] = std::string("CAL_ULD");
  m_prNames[PR_TKR_Mode] = std::string("TKR_Mode");
  m_prNames[PR_TKR_Timing] = std::string("TKR_Timing");
  m_prNames[PR_TKR_Strips] = std::string("TKR_Strips");
  m_prNames[PR_TKR_Thresh] = std::string("TKR_Thresh");
  m_prNames[PR_GNL_Mode] = std::string("GNL_Mode");
  m_prNames[PR_GNL_Timing] = std::string("GNL_Timing");
  m_prNames[PR_TRG_ROI] = std::string("TRG_ROI");
  m_prNames[PR_TRG_GEM] = std::string("TRG_GEM");
  m_prNames[PR_ACD_LCI] = std::string("ACD_LCI");
  m_prNames[PR_CAL_LCI] = std::string("CAL_LCI");
  m_prNames[PR_TKR_LCI] = std::string("TKR_LCI");

  std::vector<std::string> dbPrecincts;
  m_q->getPrecincts(dbPrecincts);
  std::vector<std::string> pclasses;
  pclasses.reserve(20);   // plenty for what we need
  int pclassIx = 0;
  for (unsigned ix = 0; ix < PR_count; ix++) {
    if (find(dbPrecincts.begin(), dbPrecincts.end(), m_prNames[ix])
        == dbPrecincts.end() ) {
      (*m_log) << MSG::ERROR << "No such precinct  (" << m_prNames[ix]
               << ")" << endreq;
      return StatusCode::FAILURE;
    }
    // If this is a latc precinct, get parm classes and store in map
    else if (( ix >= PR_firstLatc) && (ix <= PR_lastLatc)) {
      pclasses.clear();
      m_q->getParmClasses(pclasses, m_prNames[ix]);
      for (unsigned jx = 0; jx < pclasses.size(); jx++) {
        m_latcParmMap->insert(hashpair(pclasses[jx].c_str(), (int) pclassIx) );
        pclassIx++;
      }
    }
  }
  return StatusCode::SUCCESS;
}

StatusCode MootSvc::queryInterface(const InterfaceID& riid, 
                                   void** ppvInterface)
{
  if ( IID_IMootSvc == riid )  {
    // With the highest priority return the specific interface of this service
    *ppvInterface = (IMootSvc*)this;
  } else  {
    // Interface is not directly available: try out a base class
    return ConversionSvc::queryInterface(riid, ppvInterface);
  }
  addRef();
  return StatusCode::SUCCESS;
}

StatusCode MootSvc::makeMootNodes(const std::string& parent) {

  // Set up some local data structures 
  StatusCode sc = getPrecincts();
  
  std::string path(parent);
  path += "/Moot";
    
  CalibCLIDNode* node = new CalibCLIDNode(CLID_Moot);
  sc = m_provider->registerObject(path, node);
  // For now there is just one interesting node to make:
  // Make child node with additional path /Moot; then a child node
  // of that to hold a MootParmCol. These will be LATC type parms only 
  // Someday add other nodes for parm collections of other kinds
  path += "/LatcParm";
  MootParmCol* mootParmCol = new MootParmCol;

  std::string args[] = {std::string(""), std::string("")};
  unsigned long iargs[]={MOOTTYPE_parm, MOOTSUBTYPE_latcParm, 
                         (unsigned long) &m_master};
  IOpaqueAddress* pAddress;

  // Maybe this is good enough since we set the address creator in initialize?
  createAddress(MOOT_StorageType, CLID_MootParm, args, iargs,
                               pAddress);

  sc = m_dataManager->registerAddress(path, pAddress);

  sc = m_provider->registerObject(path, mootParmCol);

  if (sc.isSuccess()) m_latcPath = path;

  return sc;
}


int MootSvc::latcParmIx(const std::string& parmClass) {
  return (*m_latcParmMap)[parmClass.c_str()];
}




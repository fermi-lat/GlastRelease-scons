// $Header$

// Include files
#include "CalibDataSvc.h"
#include "CalibCLIDNode.h"
#include "CalibData/CalibTime.h"
// #include "CalibData/CalibModel.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IValidity.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/TimePoint.h"

#include "CalibData/CalibModelSvc.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static SvcFactory<CalibDataSvc> s_factory;
const ISvcFactory& CalibDataSvcFactory = s_factory;

/// Standard Constructor
CalibDataSvc::CalibDataSvc(const std::string& name,ISvcLocator* svc) :
  DataSvc(name,svc)   {
  // might also support alternative for no-network case
  declareProperty("CalibStorageType",  
                  m_calibStorageType = MYSQL_StorageType );

  // declare a property which is a list of known calibrations.
  // Have default list in one of the standard options files.  
  // User can add others.
  declareProperty("CalibNameList", m_calibList);
  declareProperty("CalibFlavorList", m_flavorList);

  declareProperty("CalibRootName",   m_calibRootName  = "Calib" ); 

  // m_rootName and m_rootCLID are declared in base class DataSvc
  m_rootName = "/" + m_calibRootName;
  m_rootCLID = CLID_DataObject;  

  m_eventTimeDefined = false;
  m_eventTime = 0;
  declareProperty("CalibInstrumentName", m_instrumentName = "LAT" );
  m_instrumentDefined = true;
}

/// Standard Destructor
CalibDataSvc::~CalibDataSvc()  {
  setDataLoader(0);
  clearStore();
}

// Service initialization
StatusCode CalibDataSvc::initialize()   {

  StatusCode sc;

  // Call base class initialisation
  sc  = DataSvc::initialize();
  if (sc.isFailure() )  return sc;

  // Set up MsgSvc
  MsgStream log(msgSvc(), name());

  // Set Data Loader
  IConversionSvc* cnv_svc;
  sc = serviceLocator()->service("DetectorPersistencySvc", cnv_svc, true);
  if (sc .isFailure() ) {
    // put something here
  }

  // Maybe...
  sc = setProperties();   
  // not sure what we'll need for properties.
  // perhaps list of known calibrations?


  sc = setDataLoader(cnv_svc);
  if (sc.isFailure() ) {
    // put something here
  }

  // Initialize the calibration data transient store
  log << MSG::DEBUG << "Storage type used is: " 
      << m_calibStorageType << endreq;
  //  log << MSG::DEBUG << "Setting CalibDataSvc root node... " << endreq;

  IAddressCreator*     calibCreator = 0;

  // Use Gaudi-supplied DetectorPersistencySvc; it's essentially
  // the same as base class PersistencySvc, which is all we need
  sc = serviceLocator()->service("DetectorPersistencySvc", calibCreator);
  
  if( sc.isFailure() ) {
    log << MSG::ERROR << "Unable to locate DetectorPersistencySvc." << endreq;
    return StatusCode::FAILURE; 
  }

  
  //   Make the root for the TDDS data
  DataObject* rootObj = new DataObject();
  sc = setRoot(m_rootName, rootObj);
  if (!sc.isSuccess() ) {
    log << MSG::ERROR << "Unable to set calib data store root." << endreq;
    delete rootObj;
    return sc;
  }

  // Create and register the next level of nodes.
  // Have one per calibration type. They are of a class trivially 
  // derived from DataObject, CalibCLIDNode.  Only additional 
  // information is CLID of child nodes.  List comes from CalibData 
  // namespace
  typedef std::vector<CalibData::CalibModelSvc::CalibPair>::const_iterator 
    PairIt;
  PairIt  pairIt;
  CalibData::CalibModelSvc svc;
  const std::vector<CalibData::CalibModelSvc::CalibPair>& pairs = 
    svc.getPairs();
  for (pairIt = pairs.begin(); pairIt != pairs.end(); pairIt++ ) {
    
    CalibCLIDNode* node = new CalibCLIDNode(pairIt->second);

    std::string calibTypePath(pairIt->first);
    sc = registerObject(calibTypePath, node);

    // Use address creator service obtained above

    IOpaqueAddress* pAddress;

    // Still have to figure out what to do about args, iargs
    unsigned long iargs[]={0, 0};

    // Set up nodes for each calibration type, default flavor
    // Create and register addresses suitable for the metadata
    // conversion service.  Ultimately, in order to find the "right"
    // set of constants,  it needs to know
    //    Calibration type, e.g. CAL Electronic gain
    //    Flavor            e.g. vanilla 
    //    Event time        validity period of constants must include this time
    //    Instrument        LAT, EM, etc.
    // We save the first two, or equivalent information, in the first
    // string parameter of a generic address
    // Consumers can use utilities in CalibData::CalibModelSvc to
    // extract fields they need
    // Event time and Instrument will be discovered by conversion service
    // when constants are requested by invoking our (CalibDataSvc) time
    // and instrument name services, resp.

    // Always do vanilla
    std::string fullpath = calibTypePath + "/vanilla";
    std::string args[] = {fullpath};
    sc = calibCreator->createAddress(m_calibStorageType, 
                                     pairIt->second,   // class id
                                     args, iargs, pAddress); 
    if (!sc.isSuccess()) {
      log << MSG::ERROR << "Unable to create Calib address" << endreq;
      }

    // A node of a specific flavor is a child of the per-calibration type
    // node for which an object was registered above.
    sc = registerAddress(fullpath, pAddress);
    if (!sc.isSuccess()) {
      log << MSG::ERROR << "Unable to register Calib address" << endreq;
    }

    // Now do the same for any requested flavors
    unsigned int ix;
    for (ix = 0; ix < m_flavorList.size(); ix++) {
      fullpath = calibTypePath + "/" + m_flavorList[ix];
      args[0] = fullpath;

      sc = calibCreator->createAddress(m_calibStorageType, 
                                       pairIt->second,
                                       args, iargs, pAddress); 
      if (!sc.isSuccess()) {
        log << MSG::ERROR << "Unable to create Calib address" << endreq;
      }

      sc = registerAddress(fullpath, pAddress);
      if (!sc.isSuccess()) {
        log << MSG::ERROR << "Unable to register Calib address" << endreq;
      }
    }    // end flavor loop 
  }      // end calibType loop

  return StatusCode::SUCCESS;
}

/// Finalize the service.
StatusCode CalibDataSvc::finalize()
{
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG << "Finalizing" << endreq;

  // Delete the associated event time
  if( 0 != m_eventTime ) delete m_eventTime; 
  m_eventTimeDefined = false;

  // Finalize the base class
  return DataSvc::finalize();
}

StatusCode CalibDataSvc::queryInterface(const IID& riid, 
				      void** ppvInterface)
{
  // With the highest priority return the specific interfaces
  // If interfaces are not directly available, try out a base class
  if ( IID_IDetDataSvc.versionMatch(riid) ) {
    *ppvInterface = (IDetDataSvc*)this;
  } else if (IID_IInstrumentName.versionMatch(riid) ) {
    *ppvInterface = (IInstrumentName*) this;
  } else if ( IID_IIncidentListener.versionMatch(riid) ) {
    *ppvInterface = (IIncidentListener*)this;
  } else {
    return DataSvc::queryInterface(riid, ppvInterface);
  }
  addRef();
  return StatusCode::SUCCESS;
}

/// Remove all data objects in the data store.
StatusCode CalibDataSvc::clearStore()   {

  MsgStream log(msgSvc(), name());
  DataSvc::clearStore();

  return StatusCode::SUCCESS;
}


/// Set the new event time 
void CalibDataSvc::setEventTime(const ITime& time) {
  m_eventTimeDefined = true;
  if (0 != m_eventTime ) delete m_eventTime; 
  m_eventTime = new CalibData::CalibTime(time);   
  MsgStream log(msgSvc(), name() );
  log << MSG::DEBUG 
      << "Event Time set to " << eventTime().hours() << endreq;
}

/// Check if the event time has been set
const bool CalibDataSvc::validEventTime() const { 
  return m_eventTimeDefined; 
}

/// Get the event time  
const ITime& CalibDataSvc::eventTime ( ) const { 
  return *m_eventTime; 
}

/// Inform that a new incident has occured
void CalibDataSvc::handle ( const Incident& inc ) { 
  MsgStream log( msgSvc(), name() );
  log << MSG::DEBUG << "New incident received" << endreq;
  log << MSG::DEBUG << "Incident source: " << inc.source() << endreq;
  log << MSG::DEBUG << "Incident type: " << inc.type() << endreq;
  return; 
}

// IInstrumentName interface
const bool CalibDataSvc::validInstrumentName() const {
  return (m_instrumentName.size() != 0); 
}

const std::string& CalibDataSvc::getInstrumentName() const {
  return m_instrumentName;
}

void CalibDataSvc::setInstrumentName(const std::string& name) {
  m_instrumentName = name;
}

StatusCode CalibDataSvc::updateObject( DataObject* toUpdate ) {

  MsgStream log( msgSvc(), name() );
  //  log << MSG::DEBUG << "Method updateObject starting" << endreq;

  // Check that object to update exists
  if ( 0 == toUpdate ) { 
    log << MSG::ERROR
	<< "There is no DataObject to update" << endreq;
    return INVALID_OBJECT; 
  }

  // Retrieve IValidity interface of object to update
  IValidity* condition = dynamic_cast<IValidity*>( toUpdate );
  if ( 0 == condition ) {
    log << MSG::WARNING
	<< "Cannot update DataObject: DataObject does not implement IValidity"
	<< endreq;
    return StatusCode::SUCCESS;
  }

  // Check that the event time has been defined
  if ( !validEventTime() ) {
    log << MSG::WARNING
	<< "Cannot update DataObject: event time undefined"
	<< endreq; 
    return StatusCode::SUCCESS;
  }

  // No need to update if condition is valid
  if ( condition->isValid( eventTime() ) ) {
    log << MSG::DEBUG 
	<< "DataObject is valid: no need to update" << endreq;
    return StatusCode::SUCCESS;
  } else {
    log << MSG::DEBUG 
	<< "DataObject is invalid: update it" << endreq;
  }

  // Now delegate update to the conversion service by calling the base class
  //  log << MSG::DEBUG 
  //      << "Delegate update to relevant conversion service" << endreq;
  StatusCode status = DataSvc::updateObject(toUpdate);
  if ( !status.isSuccess() ) {
    log << MSG::ERROR 
	<< "Could not update DataObject" << endreq; 
    if ( status == NO_DATA_LOADER )
      log << MSG::ERROR << "There is no data loader" << endreq; 
    return status;
  } 

  // Now cross-check that the new condition is valid
  condition = dynamic_cast<IValidity*>(toUpdate);
  if ( 0 == condition ) {
    log << MSG::ERROR
	<< "Updated DataObject does not implement IValidity" << endreq;
    return StatusCode::FAILURE;
  }
  if ( !condition->isValid( eventTime() ) ) {
    log << MSG::ERROR
	<< "Updated DataObject is not valid" << endreq;
    log << MSG::ERROR
	<< "Are you sure the conversion service has updated it?" << endreq;
    return StatusCode::FAILURE;
  } 

  // DataObject was successfully updated
  //  log << MSG::DEBUG << "Method updateObject exiting successfully" << endreq;
  return StatusCode::SUCCESS;
}

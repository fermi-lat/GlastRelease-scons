// $Header$

// Include files
#include "CalibDataSvc.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IValidity.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/TimePoint.h"




// Need the following to discover paths and class ids of calibration
// types.  However now CalibDataSvc depends on CalibData.  Might
// cause some problems later.
// Consider ways to
//  * localize association  of path name with class ID, perhaps in
//    CalibData (rather than CalibDataSvc)
//  * Get CalibData to do the registerObjects, or at least nicely
//    supply a vector of pairs (pathname, classID) to CalibDataSvc
//    so it can do the registerObject(s) without having to know anything.
#include "CalibData/CalibModel.h"

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

  declareProperty("CalibRootName",   m_calibRootName  = "Calib" ); 

  m_rootName = "/Calib";

  //  m_rootCLID = CLID_Catalog;
  // Probably don't need this root to do anything.  It can just be 
  // a garden variety data object.
  m_rootCLID = CLID_DataObject;  
  m_eventTimeDefined = false;
  m_eventTime = 0;
  m_instrumentName = std::string("");
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
  if (sc .isFailure() ) {
    // put something here
  }

  // Initialize the calibration data transient store
  log << MSG::DEBUG << "Storage type used is: " 
      << m_calibStorageType << endreq;
  log << MSG::DEBUG << "Setting CalibDataSvc root node... " << endreq;

  IAddressCreator*     calibCreator = 0;

  // Use Gaudi-supplied DetectorPersistencySvc; it's essentially
  // the same as base class PersistencySvc, which is all we need
  sc = serviceLocator()->service("DetectorPersistencySvc", calibCreator);
  
  if( sc.isFailure() ) {
    log << MSG::ERROR << "Unable to locate DetectorPersistencySvc." << endreq;
    return StatusCode::FAILURE; 
  }

  
  //   Make the root for the TDDS data
  DataObject rootObj = new DataObject();
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
  CalibData::PairIt  pairIt;
  for (pairIt = CalibData::pairs.begin(); pairIt++; 
       pairIt != CalibData::pairs.end(); ) {
    
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
    // two string parameters of a generic address
    // Event time and Instrument will be discovered by conversion service
    // when constants are requested by invoking our (CalibDataSvc) time
    // and instrument name services, resp.

    const std::string args[] = {calibTypePath, "/vanilla"};


    sc = calibCreator->createAddress(m_calibStorageType, 
                                     pairIt->second,   // class id
                                     args,
                                     iargs,
                                     pAddress);  //result stored here
    if (!sc.issuccess()) {
      log << MSG::ERROR << "Unable to create Calib address" 
          << endreq;
    }

    // A node of a specific flavor is a child of the per-calibration type
    // node for which an object was registered above.
    calibTypePath += "/vanilla";
    sc = registerAddress(calibTypePath, pAddress);
    if (!sc.issuccess()) {
      log << MSG::ERROR << "Unable to register Calib address" 
          << endreq;
    }
  }

  // Somebody maybe also needs to register addresses for all the
  // expected flavors of calibrations -- includes at least /vanilla
  // for each type.  Could maybe be done by CalibCnvSvc.  If in
  // CalibCnvSvc::initialize(), how to insure that it runs after
  // this set-up has been done?  
  // Might need to move above code to something that CalibCnvSvc
  // can call from its initialize() to force proper sequence.

  // OR, can ask persistency service to create the addresses (should
  // delegate to CalibCnvSvc.  This is the approach attempted above.
  //
  // what about additional flavors of standard calibration types 
  // specified by job options?
  }

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
void CalibDataSvc::setEventTime ( const ITime& time ) {
  m_eventTimeDefined = true;
  if( 0 != m_eventTime ) delete m_eventTime; 
  m_eventTime = new TimePoint( time );   
  MsgStream log( msgSvc(), name() );
  log << MSG::DEBUG 
      << "Event Time set to " << eventTime().absoluteTime() << endreq;
}

/// Check if the event time has been set
const bool CalibDataSvc::validEventTime ( ) const { 
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

void setInstrumentName(const std::string& name) {
  m_InstrumentName = name;
}


/// Update object
/// TODO: update also its ancestors in the data store if necessary
StatusCode CalibDataSvc::updateObject( DataObject* toUpdate ) {

  MsgStream log( msgSvc(), name() );
  log << MSG::DEBUG << "Method updateObject starting" << endreq;

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

  // TODO: before loading updated object, update HERE its parent in data store

  // Now delegate update to the conversion service by calling the base class
  log << MSG::DEBUG 
      << "Delegate update to relevant conversion service" << endreq;
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
  log << MSG::DEBUG << "Method updateObject exiting successfully" << endreq;
  return StatusCode::SUCCESS;

}

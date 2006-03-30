// $Header$

// Include files
#include "CalibDataSvc.h"
#include "CalibCLIDNode.h"
#include "CalibData/CalibTime.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IValidity.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/TimePoint.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/TopLevel/Event.h"
#include "astro/JulianDate.h"

#include "CalibData/CalibModelSvc.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static SvcFactory<CalibDataSvc> s_factory;
const ISvcFactory& CalibDataSvcFactory = s_factory;

/// Standard Constructor
CalibDataSvc::CalibDataSvc(const std::string& name,ISvcLocator* svc) :
  DataSvc(name,svc), m_useEventTime(true), m_timeSourceEnum(TIMESOURCEnone) {
  // might also support alternative for no-network case
  declareProperty("CalibStorageType",  
                  m_calibStorageType = MYSQL_StorageType );

  // declare a property which is a list of known calibrations.
  // Have default list in one of the standard options files.  
  // User can add others.
  declareProperty("CalibNameList", m_calibList);
  declareProperty("CalibFlavorList", m_flavorList);
  declareProperty("CalibRootName",   m_calibRootName  = "Calib" ); 
  declareProperty("UseEventTime", m_useEventTime = true);

  // m_rootName and m_rootCLID are declared in base class DataSvc
  m_rootName = "/" + m_calibRootName;
  m_rootCLID = CLID_DataObject;  

  m_eventTimeDefined = false;
  m_eventTime = 0;
  m_nEvent = 0;
  declareProperty("CalibInstrumentName", m_instrumentName = "LAT" );

  // choices could be "data", "clock", "mc", "none", "digi"
  //  declareProperty("CalibTimeSource", m_timeSource = "none" );
  declareProperty("CalibTimeSource", m_timeSource = "data" );
  m_instrumentDefined = true;

  // Still more properties.  Used only by fake clock, if active
  /*  declareProperty( "startTime",  
      m_startTimeAsc = "2003-1-10_00:20");   */
  declareProperty( "startTime",  
                   m_startTimeAsc = "2003-1-10 00:20");

  declareProperty( "delayTime",  m_delayTime = 2000);

}

/// Standard Destructor
CalibDataSvc::~CalibDataSvc()  {
  setDataLoader(0);
  clearStore();
}

// Service initialization
StatusCode CalibDataSvc::initialize()   {

  StatusCode sc;

  sc  = DataSvc::initialize();     // Call base class initialisation
  if (sc.isFailure() )  return sc;

  // Set up MsgSvc, Data Loader
  MsgStream log(msgSvc(), name());
  IConversionSvc* cnv_svc;
  sc = serviceLocator()->service("DetectorPersistencySvc", cnv_svc, true);
  if (sc .isFailure() ) {
    log << MSG::ERROR << "Unable to find DetectorPersistencySvc " << endreq;
    return sc;
  }

  // Need event data service for timestamp stuff
  sc = serviceLocator()->service("EventDataSvc", m_eventSvc, true);
  if (sc .isFailure() ) {
    log << MSG::ERROR << "Unable to find EventDataSvc " << endreq;
    return sc;
  }


  sc = setDataLoader(cnv_svc);
  if (sc.isFailure() ) {
    log << MSG::ERROR << "Unable to set data loader " << endreq;
    return sc;
  }
  sc = setProperties();   

  // Initialize the calibration data transient store
  log << MSG::DEBUG << "Storage type used is: " 
      << m_calibStorageType << endreq;

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

  // Get ready to listen for BeginEvent
  if (m_timeSource != "none") {
    if (m_timeSource == "data") {
      m_timeSourceEnum = TIMESOURCEdata;
    }
    else if (m_timeSource == "mc") {
      m_timeSourceEnum = TIMESOURCEmc;
    }
    else if (m_timeSource == "digi") {
      m_timeSourceEnum = TIMESOURCEdigi;
    }
    else if (m_timeSource == "clock") {
      m_timeSourceEnum = TIMESOURCEclock;

      //  set start and increment parameters also
      unsigned int underpos = m_startTimeAsc.find("_");
      if (underpos < m_startTimeAsc.size()) {
        m_startTimeAsc.replace(underpos, 1, " ");
      }
      m_startTime = facilities::Timestamp(m_startTimeAsc).getClibTime();

    }
    else {
      log << MSG::WARNING << "Unsupported time source " << m_timeSource 
          << " will be set to  'none' " << endreq;
      m_timeSource = std::string("none");
    }

    IIncidentSvc* incSvc;
    StatusCode sc = service("IncidentSvc", incSvc, true);
    if (sc.isSuccess() ) {
      int priority = 100;
      incSvc->addListener(this, "BeginEvent", priority);
    }
    else {
      log << MSG::ERROR << "Unable to find IncidentSvc" << endreq;
      return sc;
    }
  }
  // Make flavor nodes in the calibration TDS
  return makeFlavorNodes(calibCreator, &log);
}

// Create and register the next level of nodes.
// Have one per calibration type. They are of a class trivially 
// derived from DataObject, CalibCLIDNode.  Only additional 
// information is CLID of child nodes.  List comes from CalibData 
// namespace
StatusCode CalibDataSvc::makeFlavorNodes(IAddressCreator*  calibCreator,
                                         MsgStream* log) {
  typedef std::vector<CalibData::CalibModelSvc::CalibPair>::const_iterator 
    PairIt;
  PairIt  pairIt;
  StatusCode sc;
  CalibData::CalibModelSvc svc;
  const std::vector<CalibData::CalibModelSvc::CalibPair>& pairs = 
    svc.getPairs();
  for (pairIt = pairs.begin(); pairIt != pairs.end(); pairIt++ ) {
    
    CalibCLIDNode* node = new CalibCLIDNode(pairIt->second);

    std::string calibTypePath(pairIt->first);
    sc = registerObject(calibTypePath, node);

    // Still have to figure out what to do about args, iargs
    unsigned long iargs[]={0, 0};
    IOpaqueAddress* pAddress;

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
      (*log) << MSG::ERROR 
          << "Unable to create Calib address with path " << fullpath << endreq;
    }

    // A node of a specific flavor is a child of the per-calibration type
    // node for which an object was registered above.
    sc = registerAddress(fullpath, pAddress);
    if (!sc.isSuccess()) {
      (*log) << MSG::ERROR << "Unable to register Calib address with path" 
          << fullpath << endreq;
    }
    // Now do the same for any requested flavors
    unsigned int ix;

    for (ix = 0; ix < m_flavorList.size(); ix++) {
      // Don't redo vanilla, if present in user list
      if (m_flavorList[ix] == std::string("vanilla")) continue;

      fullpath = calibTypePath + "/" + m_flavorList[ix];
      args[0] = fullpath;

      sc = calibCreator->createAddress(m_calibStorageType, 
                                       pairIt->second, args, iargs, pAddress); 
      if (!sc.isSuccess()) {
        (*log) << MSG::ERROR << "Unable to create Calib address with path " 
            << fullpath << endreq;
      }
      sc = registerAddress(fullpath, pAddress);
      if (!sc.isSuccess()) {
        (*log) << MSG::ERROR << "Unable to register Calib address with path " 
            << fullpath << endreq;
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

StatusCode CalibDataSvc::queryInterface(const InterfaceID& riid, void** ppvInterface)
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
  //  MsgStream log(msgSvc(), name());

  DataSvc::clearStore();
  return StatusCode::SUCCESS;
}


/// Set the new event time 
void CalibDataSvc::setEventTime(const ITime& time) {
  m_eventTimeDefined = true;
  if (0 != m_eventTime ) delete m_eventTime; 
  CalibData::CalibTime *pTime = new CalibData::CalibTime(time);   
  m_eventTime = pTime;

  if ((m_nEvent < 100) || (m_nEvent == ((m_nEvent/100) * 100) ) ) {
    MsgStream log(msgSvc(), name() );
    log << MSG::DEBUG 
        << "Event Time set to " 
        << pTime->getString() << endreq;
  }

}

/// Check if the event time has been set
const bool CalibDataSvc::validEventTime() const { 
  return (m_eventTimeDefined); 
}

/// Get the event time  
const ITime& CalibDataSvc::eventTime ( )  const { 
  using CalibData::CalibTime;

  static CalibTime   badTime(facilities::Timestamp(0,0));
  if (m_eventTime) {
    return *m_eventTime;
  }
  else return badTime;
}

/// Inform that a new incident has occured
void CalibDataSvc::handle ( const Incident& inc ) { 

  if ((inc.type() == "BeginEvent") && 
      ((m_timeSourceEnum == TIMESOURCEdata) ||
       (m_timeSourceEnum == TIMESOURCEmc) ||
       (m_timeSourceEnum == TIMESOURCEdigi) ||
       (m_timeSourceEnum == TIMESOURCEclock) ) ) {
    if ((m_nEvent < 100) || (m_nEvent == ((m_nEvent/100) * 100) ) ) {
      MsgStream log( msgSvc(), name() );
      log << MSG::DEBUG << "New incident received" << endreq;
      log << MSG::DEBUG << "Incident source: " << inc.source() << endreq;
      log << MSG::DEBUG << "Incident type: " << inc.type() << endreq;
    }
    
    m_nEvent++;
    m_newEvent = true;

  }
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
  // Check that object to update exists
  if ( 0 == toUpdate ) { 
    log << MSG::ERROR
	<< "There is no DataObject to update" << endreq;
    return INVALID_OBJECT; 
  }

  // Update timestamp if necessary
  updateTime();

  // Retrieve IValidity interface of object to update
  IValidity* condition = dynamic_cast<IValidity*>( toUpdate );
  if ( 0 == condition ) {
    log << MSG::WARNING
	<< "Cannot update DataObject: DataObject does not implement IValidity"
	<< endreq;
    return StatusCode::SUCCESS;
  }

  if (m_useEventTime) {
    // Check that the event time has been defined
    if ( !validEventTime() ) {
      log << MSG::WARNING
          << "Cannot update DataObject: event time undefined"
          << endreq; 
      return StatusCode::SUCCESS;
    }

    // No need to update if condition is valid
    if ( condition->isValid( eventTime() ) ) {
      if ((m_nEvent < 100) || (m_nEvent == ((m_nEvent/100) * 100) ) ) {
        log << MSG::DEBUG 
            << "DataObject is valid: no need to update" << endreq;
      }
      return StatusCode::SUCCESS;
    } else {
      log << MSG::DEBUG 
          << "DataObject is invalid: update it" << endreq;
    }
  }
  // *** May want to check for analogous condition for not-using-event-time
  // case as well

  // Now delegate update to the conversion service by calling the base class
  StatusCode status = DataSvc::updateObject(toUpdate);
  if ( !status.isSuccess() ) {
    log << MSG::ERROR 
	<< "Could not update DataObject" << endreq; 
    if ( status == NO_DATA_LOADER )
      log << MSG::ERROR << "There is no data loader" << endreq; 
    return status;
  } 

  // If we're not using event time, nothing to check
  if (!m_useEventTime)  return StatusCode::SUCCESS;

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
  return StatusCode::SUCCESS;
}

StatusCode CalibDataSvc::loadObject(IConversionSvc* pLoader, 
                                    IRegistry* pRegistry) {
  if (m_newEvent) {
    updateTime();
    m_newEvent = false;
  }
  return DataSvc::loadObject(pLoader, pRegistry);

}

// For timeSource = "data", "mc", "digi" or "clock"
StatusCode  CalibDataSvc::updateTime() {
  using CalibData::CalibTime;

  //  MsgStream log( msgSvc(), name() );


  if (!m_useEventTime) return StatusCode::SUCCESS;
  if (m_timeSource == "none") return StatusCode::SUCCESS;

  StatusCode sc = StatusCode::FAILURE;

  if (m_newEvent) {
    // Fetch the time using requested fetch mechanism
    switch (m_timeSourceEnum) {
    case TIMESOURCEdata:
    case TIMESOURCEdigi: 
      {

      sc = fetchEventTime();

      // no longer need to handle separately with these calls:
      // sc = fetchDataTime();
      // sc = fetchMcTime();   ???
      // sc = fetchDigiTime();

      break;
    }
    case TIMESOURCEmc: {
      sc = fetchMcTime();
      //      sc = fetchEventTime();
      break;
    }
    case TIMESOURCEclock: {
      sc = fetchFakeClockTime();
      break;
    }
    default: {
      sc = StatusCode::FAILURE;
    }

    }             // end of switch

    if (sc.isSuccess()) setEventTime(m_time);

    m_newEvent = false;
    return sc;
  }
  return StatusCode::SUCCESS;
}
/* ---    */

StatusCode CalibDataSvc::fetchEventTime() {

  MsgStream log(msgSvc(), name());

  static const facilities::Timestamp missionStart("2001-1-1 00:00");
  static const unsigned missionSec = (unsigned) missionStart.getClibTime();
  static const int missionNano = missionStart.getNano();

  SmartDataPtr<Event::EventHeader> eventHeader(m_eventSvc, "/Event");
  
  if (!eventHeader) {
    log << MSG::ERROR << "Unable to retrieve event timestamp for digis" 
        << endreq;
    return StatusCode::FAILURE;
  }
  double fromMissionStart = (eventHeader->time()).time();
  unsigned fromSec = (unsigned) fromMissionStart;
  unsigned fromNano = (unsigned) ((fromMissionStart - fromSec) * 1000000000);

  if ((m_nEvent < 100) || (m_nEvent == ((m_nEvent/100) * 100) ) ) {

    log << MSG::DEBUG << "event time in seconds/nano from mission start: " 
        << fromSec 
        << "/" << fromNano << endreq;
  }
  facilities::Timestamp absTime(missionSec + fromSec, missionNano + fromNano);
  m_time = CalibData::CalibTime(absTime);

  return StatusCode::SUCCESS; 
}


StatusCode CalibDataSvc::fetchMcTime() {

  MsgStream log(msgSvc(), name());

  //  static const facilities::Timestamp missionStart("2001-1-1 00:00");
  static const astro::JulianDate 
    julianStart = astro::JulianDate::missionStart();
  static const facilities::Timestamp missionStart(julianStart);

  static const unsigned missionSec = (unsigned) missionStart.getClibTime();
  static const int missionNano = missionStart.getNano();

  SmartDataPtr<Event::MCEvent> mcHeader(m_eventSvc, "/Event/MC");
  if (!mcHeader) {
    log << MSG::ERROR << "Unable to retrieve mc event timestamp " << endreq;

    return StatusCode::FAILURE;
  }
  double fromMissionStart = (mcHeader->time()).time();
  unsigned fromSec = (unsigned) fromMissionStart;
  unsigned fromNano = (unsigned) ((fromMissionStart - fromSec) * 1000000000);

  if ((m_nEvent < 100) || (m_nEvent == ((m_nEvent/100) * 100) ) ) {
    log << MSG::DEBUG << "(mc) seconds/nano from mission start: " << fromSec 
        << "/" << fromNano << endreq;
  }

  facilities::Timestamp absTime(missionSec + fromSec, missionNano + fromNano);
  m_time = CalibData::CalibTime(absTime);

  return StatusCode::SUCCESS; 
}


StatusCode CalibDataSvc::fetchFakeClockTime() {

  static unsigned int eventNumber = 0;
  static long seconds = 0;
  static long nano    = 0;

  static long seconds_incr = m_delayTime / 1000;
  static long nano_incr = (m_delayTime - 1000*seconds_incr) * 1000000;

  seconds += seconds_incr;
  nano += nano_incr;
  while (nano >= 1000000000) {
    nano -= 1000000000;
    seconds++;
  }
  // Set the event time.  m_delay has been specified in milliseconds
    //  double incr = eventNumber*m_delayTime;
    //  long   seconds_incr = (long) (incr / 1000);
    //  long   nano_incr = (long) ((incr - (1000*seconds_incr)) * 1000000);

  std::string timestring;
  try {
    facilities::Timestamp time(m_startTime + seconds, nano);
    timestring = time.getString();
    m_time = CalibData::CalibTime(time);
  }
  catch (facilities::BadTimeInput  ex) {
    std::cerr << "Caught facilities::BadTimeInput exception with complaint "
              << ex.complaint << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::cerr.flush();
    exit(1);
  }

  if ((m_nEvent < 100) || (m_nEvent == ((m_nEvent/100) * 100) ) ) {
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << "Fake Clock Event number: " 
        << eventNumber << endreq;
    log << MSG::DEBUG << "Event time: "
        << timestring
        << endreq; 
    log << MSG::DEBUG << "Event time (hours) " << m_time.hours() << endreq;
  }
  eventNumber++;

  return StatusCode::SUCCESS; 
}

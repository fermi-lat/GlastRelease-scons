//$Header$
#include <string>
#include <cstdio>
#include <stdexcept>

#include "CalibMySQLCnvSvc.h"
#include "CalibSvc/ICalibRootSvc.h"   // for def. of CALIBROOT_StorageType
#include "calibUtil/Metadata.h"

#include "CalibData/CalibBase.h"
#include "CalibData/CalibTime.h"
#include "CalibData/CalibModelSvc.h"

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/GenericAddress.h"
#include "GaudiKernel/IConverter.h"
#include "GaudiKernel/IDetDataSvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IValidity.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"

/// Instantiation of a static factory to create instances of this service
//static SvcFactory<CalibMySQLCnvSvc>          CalibMySQLCnvSvc_factory;
//const ISvcFactory& CalibMySQLCnvSvcFactory = CalibMySQLCnvSvc_factory;
DECLARE_SERVICE_FACTORY(CalibMySQLCnvSvc);

// Local utility to translate calibration quality list to bit map
namespace {
  unsigned int toQualityMask(std::vector<std::string>& qualities) {
    using calibUtil::Metadata;

    unsigned int mask = 0;
    unsigned n = qualities.size();

    for (unsigned i = 0; i < n; i++) {
      std::string iString = qualities[i];
      if (iString.size() < 3) continue;
      iString.resize(3);
      if (iString == "PRO") mask |= Metadata::LEVELProd;
      else if (iString =="DEV") mask |= Metadata::LEVELDev;
      else if (iString =="TES") mask |= Metadata::LEVELTest;
      else if (iString =="SUP") mask |= Metadata::LEVELSuperseded;
    }
    return mask;
  }
}

CalibMySQLCnvSvc::CalibMySQLCnvSvc( const std::string& name, ISvcLocator* svc)
  : ConversionSvc (name, svc, MYSQL_StorageType)
    , m_meta(0), m_useEventTime(true),m_enterTimeStart(0), m_enterTimeEnd(0),
    m_qualityMask(0)
{
  declareProperty("Host", m_host = "*");
  declareProperty("UseEventTime", m_useEventTime = true);
  declareProperty("EnterTimeEnd", m_enterTimeEndString = std::string("") );
  declareProperty("EnterTimeStart", m_enterTimeStartString = std::string("") );
  declareProperty("DbName", m_dbName = std::string("calib") );
  declareProperty("QualityList", m_qualityList);
  declareProperty("CrashOnError", m_crash = true);
}

CalibMySQLCnvSvc::~CalibMySQLCnvSvc(){ }

StatusCode CalibMySQLCnvSvc::initialize()
{
  // Initialize base class
  StatusCode sc = ConversionSvc::initialize();
  if ( !sc.isSuccess() ) return sc;

  MsgStream log(msgSvc(), "CalibMySQLCnvSvc" );
  log << MSG::INFO << "Specific initialization starting" << endreq;

  // Locate the Calib Data Service.  Since it inherits from DataSvc
  // it has to implement IDataProviderSvc
  IDataProviderSvc* pCDS = 0;
  sc = serviceLocator()->getService 
      ("CalibDataSvc",  IDataProviderSvc::interfaceID(), (IInterface*&)pCDS);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR << "Could not locate CalibDataSvc" << endreq;
    return sc;
  }

  // Set the CalibDataSvc as data provider service
  sc = setDataProvider(pCDS);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR << "Could not set data provider" << endreq;
    return sc;
  }

  // Query the IDetDataSvc and IInstrumentName interfaces of the 
  // calib data service
  sc = pCDS->queryInterface(IDetDataSvc::interfaceID(), (void**) &m_detDataSvc);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Cannot query IDetDataSvc interface of CalibDataSvc" << endreq;
    return sc;
  } else {
    log << MSG::DEBUG << "Retrieved IDetDataSvc interface of CalibDataSvc" 
	<< endreq;
  }

  sc = pCDS->queryInterface(IID_IInstrumentName, (void**) &m_instrSvc);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Cannot query IInstrumentName interface of CalibDataSvc" << endreq;
    return sc;
  } else {
    log << MSG::DEBUG 
        << "Retrieved IInstrumentNameSvc interface of CalibDataSvc" << endreq;
  }

  // Locate IConversionSvc interface of the DetectorPersistencySvc
  sc = serviceLocator()->service 
    ("DetectorPersistencySvc", m_detPersSvc, true);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Cannot locate IConversionSvc interface of DetectorPersistencySvc"
	<< endreq;
    return sc;
  } else {
    log << MSG::DEBUG 
	<< "Retrieved IConversionSvc interface of DetectorPersistencySvc"
	<< endreq;
  }
  
  // Query the IAddressCreator interface of the detector persistency service
  IAddressCreator* iAddrCreator;
  sc = m_detPersSvc->queryInterface(IAddressCreator::interfaceID(), 
				    (void**) &iAddrCreator);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 
	<< "Cannot query IAddressCreator interface of DetectorPersistencySvc" 
	<< endreq;
    return sc;
  } else {
    log << MSG::DEBUG 
	<< "Retrieved IAddressCreator interface of DetectorPersistencySvc" 
	<< endreq;
  }
  //  log << MSG::DEBUG 
  //    << "Set it as the address creator of the CalibMySQLCnvSvc" << endreq;
  sc = setAddressCreator(iAddrCreator);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 	<< "Cannot set the address creator" << endreq;
    return sc;
  }

  // Get properties from the JobOptionsSvc
  sc = setProperties();
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR << "Could not set jobOptions properties" << endreq;
    return sc;
  }
  log << MSG::DEBUG << "Properties were read from jobOptions" << endreq;

  // Translate list of calibration quality names to bit mask form used 
  // by calibUtil::Metadata::findBest     Defaults to PROD + DEV for now
  // (that was old fixed value)
  m_qualityMask = toQualityMask(m_qualityList);
  if (!m_qualityMask) {
    m_qualityMask = calibUtil::Metadata::LEVELProd |
      calibUtil::Metadata::LEVELDev;
  }


  if (!m_useEventTime) {  // special diagnostic mode
    std::string defStart("2003-09-01");
    if (m_enterTimeStartString.size() == 0) {
      // set some  default
      m_enterTimeStartString = defStart;
    }
    log << MSG::INFO << "Calibrations created at or after " 
        << m_enterTimeStartString << " will be fetched" << endreq;
    try {
      m_enterTimeStart = new facilities::Timestamp(m_enterTimeStartString);
    }
    catch (facilities::BadTimeInput ex) {
      log << MSG::ERROR << "Enter time no good:  " << ex.complaint << endreq;
      log << "Using default " << defStart << endreq;
      m_enterTimeStartString = defStart;
      m_enterTimeStart = new facilities::Timestamp(m_enterTimeStartString);
    }
    if (m_enterTimeEndString.size()) {
      try {
        m_enterTimeEnd = new facilities::Timestamp(m_enterTimeEndString);
      }
      catch (facilities::BadTimeInput ex) {
        log << MSG::ERROR << "Enter end time no good: " << ex.complaint 
            << endreq;
        log << "Using default of forever " << endreq;
        m_enterTimeEnd = 0;
      }
    }
    // There is no easy way to tell data service about this, so user
    // is expected to supply yet another job options parameter to CalibDataSvc
    // to tell it not to check for valid event time.  Ugh.
  }
  // Make a calibUtil::Metadata instance 
  // Conceivably, could start up a different conversion service, depending 
  // on job options parameters, which would look very much like this one
  // except for having a different way to access metadata.
  m_meta = new calibUtil::Metadata(m_host, "*", m_dbName);

  if (!m_meta) {
    log << MSG::ERROR << "Could not open connection to metadata dbs" << endreq;
    return MSG::ERROR;
  }

  // Arrange to be woken up once per event
  IIncidentSvc* incSvc;
  sc = service("IncidentSvc", incSvc, true);
  if (sc.isSuccess() ) {
    int priority = 100;
    incSvc->addListener(this, "EndEvent", priority);
  }
  else {
    log << MSG::ERROR << "Unable to find IncidentSvc" << endreq;
    return sc;
  }


  // Probably should get this value from job options. 
  // Now we do.  See m_qualityMask, m_qualityList
  //  m_calibLevelMask = calibUtil::Metadata::LEVELProd + 
  //    calibUtil::Metadata::LEVELDev;

  log << MSG::INFO << "Specific initialization completed" << endreq;
  return sc;
}


StatusCode CalibMySQLCnvSvc::finalize()
{
  MsgStream log(msgSvc(), "CalibMySQLCnvSvc");
  log << MSG::DEBUG << "Finalizing" << endreq;
  delete m_meta;
  m_meta = 0;
  return ConversionSvc::finalize();
}


StatusCode CalibMySQLCnvSvc::queryInterface(const InterfaceID& riid, 
                                            void** ppvInterface)
{
  if ( IID_ICalibMetaCnvSvc == riid )  {
    // With the highest priority return the specific interface of this service
    *ppvInterface = (ICalibMetaCnvSvc*)this;
  } else  {
    // Interface is not directly available: try out a base class
    return ConversionSvc::queryInterface(riid, ppvInterface);
  }
  addRef();
  return StatusCode::SUCCESS;
}


/// Create a transient representation from another representation of an object.
/// Overloaded from ConversionSvc because CalibMySQLCnvSvc has no converters.
/// (The typical conversion service delegates this to an appropriate converter)
StatusCode CalibMySQLCnvSvc::createObj (IOpaqueAddress* pAddress, 
					DataObject*&    refpObject ) {

  MsgStream log(msgSvc(), "CalibMySQLCnvSvc" );
  // log << MSG::DEBUG << "Method createObj starting" << endreq;

  if (( (!m_detDataSvc->validEventTime()) && m_useEventTime ) ||
      ( !m_instrSvc->validInstrumentName())  )   {
    log << MSG::ERROR
	<< "Cannot create DataObject: event time or instrument undefined"
	<< endreq; 
    return StatusCode::FAILURE;
  } 


  // Create the object according to calib type, flavor, time, instrument, clid.
  // Notice that the CalibMySQLCnvSvc has no converters of its own:
  // object creation is delegated to another CnvSvc via a temporary address
  // The IOpaqueAddress specifies calibration type and specific flavor.
  // The secondary storage type is always discovered dynamically
  StatusCode sc;

  sc = createCalib(refpObject,
                   pAddress->par()[0],   // full path
                   m_detDataSvc->eventTime(),
                   m_instrSvc->getInstrumentName(),
                   pAddress->clID(),
                   pAddress->registry() );

  if ( !sc.isSuccess() ) {
    log << MSG::ERROR << "Could not create calib DataObject" << endreq;
  }
  log << MSG::DEBUG << "Method createObj exiting" << endreq;
  return sc;
}

/// Resolve the references of the created transient object.
/// (Actually, don't, because this operation isn't supported, nor is
/// it needed for the conversion service.)
/// Overloaded from ConversionSvc because CalibMySQLCnvSvc has no converters.
StatusCode CalibMySQLCnvSvc::fillObjRefs(IOpaqueAddress* /*pAddress*/, 
                                         DataObject*     /*pObject */ ) {
  MsgStream log(msgSvc(), "CalibMySQLCnvSvc" );
  // log << MSG::DEBUG << "Method fillObjRefs is not implemented" << endreq;
  return StatusCode::SUCCESS;
}
  

/// Update a transient representation from another representation of an object.
/// Overloaded from ConversionSvc because CalibMySQLCnvSvc has no converters.
StatusCode CalibMySQLCnvSvc::updateObj(IOpaqueAddress* pAddress, 
                                       DataObject*     pObject  ) {

  using facilities::Timestamp;

  // Don't update when we're using enter time
  if (!m_useEventTime) return StatusCode::SUCCESS;

  MsgStream log(msgSvc(), "CalibMySQLCnvSvc" );
  //  log << MSG::DEBUG << "Method updateObj starting" << endreq;

  // Check that the event time has been defined
  // The event time is obtained from the ConditionDataSvc
  if ( !m_detDataSvc->validEventTime() ) {
    log << MSG::ERROR
	<< "Cannot update DataObject: event time undefined"
	<< endreq; 
    return StatusCode::FAILURE;
  } else {
    if (m_useEventTime) {
      std::string tString  = 
        CalibData::CalibTime(m_detDataSvc->eventTime()).getString();
      log << MSG::DEBUG
          << "Event time: ";
      log <<  tString << endreq;
    }
  }
  if ( !m_instrSvc->validInstrumentName() ) {
    log << MSG::ERROR
	<< "Cannot update DataObject: instrument name undefined"
	<< endreq; 
    return StatusCode::FAILURE;
  }

  // Update the object according to calibration type, flavor, time, 
  // instrument, clid, 
  // Notice that the CalibMySQLCnvSvc has no converters of its own:
  // object creation is delegated to another CnvSvc via a temporary address
  // The IOpaqueAddress specifies calibration type and flavor
  // The secondary storage type is always discovered dynamically
  if( 0 == pObject ) {
    log << MSG::ERROR << "There is no object to update" << endreq;
    return StatusCode::FAILURE;
  }
  IValidity* pValidity = dynamic_cast<IValidity*>(pObject);
  if (0 == pValidity) {
    log << MSG::WARNING
	<< "Object to update does not implement IValidity: assume up-to-date" 
	<< endreq;
    return StatusCode::SUCCESS;
  }

  std::string tString  = 
    CalibData::CalibTime(pValidity->validSince()).getString();

  log << MSG::DEBUG << "Old calib DataObject was valid since "
      << tString << " till ";

  log << CalibData::CalibTime(pValidity->validTill()).getString() << endreq;

  StatusCode sc;
  if (pValidity->isValid(m_detDataSvc->eventTime())) { // done
    return StatusCode::SUCCESS;
  }

  sc = updateCalib(pObject, pAddress->par()[0], 
                   m_detDataSvc->eventTime(), m_instrSvc->getInstrumentName(),
                   pAddress->clID(), pAddress->registry() );
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR << "Could not update calib DataObject" << endreq;
    return sc;
  }

  // Last, check that everything is OK
  pValidity = dynamic_cast<IValidity*>(pObject);
  if ( 0 == pValidity ) {
    log << MSG::ERROR 
	<< "Updated object does not implement IValidity" << endreq;
    return StatusCode::FAILURE;
  }
  log << MSG::DEBUG << "New calib DataObject is valid (year only)  since "
  << pValidity->validSince().year(false) 
  << " till " << pValidity->validTill().year(false)  
  << endreq;

  //  log << MSG::DEBUG << "Method updateObj exiting" << endreq;
  return StatusCode::SUCCESS;  
}

/// Update the references of an updated transient object. [actually, don't.
/// Calib data doesn't have any inter-object references.]
/// Overloaded from ConversionSvc because CalibMySQLCnvSvc has no converters.
StatusCode CalibMySQLCnvSvc::updateObjRefs (IOpaqueAddress* /*pAddress*/, 
                                            DataObject*     /*pObject */ ) {
  MsgStream log(msgSvc(), "CalibMySQLCnvSvc" );
  //  log << MSG::WARNING << "Method updateObjRefs is not implemented" 
  // << endreq;
  return StatusCode::SUCCESS;
}
  

/// Convert a transient object to a requested representation. Not implemented.
/// Overloaded from ConversionSvc because CalibMySQLCnvSvc has no converters.
StatusCode CalibMySQLCnvSvc::createRep(DataObject* /*pObject*/,
				       IOpaqueAddress*& /*refpAddress*/ ) {

  MsgStream log(msgSvc(), "CalibMySQLCnvSvc" );
  // log << MSG::WARNING << "Method createRep is not implemented" << endreq;
  return StatusCode::SUCCESS;
}
  

/// Resolve the references of a converted object. [actually, don't.
/// Calib data doesn't have any inter-object references.]
/// Overloaded from ConversionSvc because CalibMySQLCnvSvc has no converters.
StatusCode CalibMySQLCnvSvc::fillRepRefs (IOpaqueAddress* /*pAddress*/, 
                                          DataObject*     /*pObject */ ) {
  MsgStream log(msgSvc(), "CalibMySQLCnvSvc" );
  //  log << MSG::WARNING << "Method fillRepRefs is not implemented" << endreq;
  return StatusCode::SUCCESS;
}
  

/// Update a converted representation of a transient object.
/// Overloaded from ConversionSvc because CalibMySQLCnvSvc has no converters.
StatusCode CalibMySQLCnvSvc::updateRep (IOpaqueAddress* /*pAddress*/, 
                                        DataObject*     /*pObject */ ) {
  MsgStream log(msgSvc(), "CalibMySQLCnvSvc" );
  // log << MSG::WARNING << "Method updateRep is not implemented" << endreq;
  return StatusCode::SUCCESS;
}
  

/// Update the references of an already converted object.
/// Overloaded from ConversionSvc because CalibMySQLCnvSvc has no converters.
/// Don't do anything because calib objects have no inter-object references.
StatusCode CalibMySQLCnvSvc::updateRepRefs (IOpaqueAddress* /*pAddress*/, 
                                            DataObject*     /*pObject */ ) {
  MsgStream log(msgSvc(), "CalibMySQLCnvSvc");
  //log << MSG::WARNING << "Method updateRepRefs is not implemented" << endreq;
  return StatusCode::SUCCESS;
}
  
/// Overload ConversionSvc implementation of createAddress.  
/// Create an address using explicit arguments to identify a single object.
/// Par[0] is full path in calibration TDS
StatusCode CalibMySQLCnvSvc::createAddress(long svc_type,
                                           const CLID& clid,
                                           const std::string* par, 
                                           const unsigned long* /*ipar*/,
                                           IOpaqueAddress*& refpAddress ) {

  // First check that requested address is of type MYSQL_StorageType
  MsgStream log(msgSvc(), "CalibMySQLCnvSvc" );
  if ( svc_type!= MYSQL_StorageType ) {
    log << MSG::ERROR 
	<< "Cannot create addresses of type " << (int)svc_type 
	<< " which is different from " << (int)MYSQL_StorageType 
	<< endreq;
    return StatusCode::FAILURE;
  }
  
  refpAddress = new GenericAddress( MYSQL_StorageType, 
				    clid, 
                                    par[0]);

  return StatusCode::SUCCESS;
}


/// Create a calib DataObject by calib type, flavor, time and instrument.
/// This method does not register DataObject in the transient data store,
/// [but may register TDS addresses for its children if needed (e.g. Catalog).
///   - what's all this about? ]
/// The string storage type is discovered at runtime in the metadata dbs
/// The entry name identifies a condition amongst the many in the string.
/// Implementation:
/// - create a temporary address containing storage type and classID;
/// - dispatch to appropriate conversion service according to storage type;
/// - this will dispatch to appropriate converter according to CLID
///   (CalibMySQLCnvSvc has no converters of its own).
StatusCode CalibMySQLCnvSvc::createCalib(DataObject*&       refpObject,
                                         const std::string& fullpath,
                                         const Gaudi::Time&       time,
                                         const std::string& instrName,
                                         const CLID&        classID,
                                         IRegistry*         entry)
{
  MsgStream log(msgSvc(), "CalibMySQLCnvSvc" );

  // Look up calib object in the Metadata database
  std::string flavor = CalibData::CalibModelSvc::getFlavor(fullpath);
  std::string cType = CalibData::CalibModelSvc::getCalibType(fullpath);

  // ..and extra special munging for test
  if (std::string("Test") == cType.substr(0, 4)) {
    cType = std::string("Test_Gen");
  }

  unsigned int ser;
  m_metaConnected = true;  // one way or another, will be making connection

  calibUtil::Metadata::eRet ret;
  if (m_useEventTime) {
    ret = m_meta->findBest(&ser, cType, CalibData::CalibTime(time),
                           m_qualityMask, instrName, flavor);
  }
  else {
    ret = m_meta->findSoonAfter(&ser, cType, m_enterTimeStart, 
                                m_enterTimeEnd, m_qualityMask, 
                                instrName, flavor);
  }
  MSG::Level msgLevel = (m_crash) ? MSG::FATAL : MSG::ERROR;
  if (ret != calibUtil::Metadata::RETOk) {
    log << msgLevel << "Could not access MySQL database" << endreq;
    //    if (m_crash) throw std::runtime_error("Could not acces MySQL dbs");
    if (m_crash) {
      log << msgLevel << std::endl << "Exiting... " << std::endl << endreq;
      exit(1);
    }
    return StatusCode::FAILURE;
  }
  else if (ser == 0) {  // no error, but no appropriate calib was found
    log << msgLevel
        << std::endl << "No appropriate calibration of (type, flavor) ("
        << cType << ", " << flavor << ")" << endreq;
    if (m_crash) {
      log << msgLevel << std::endl << "Exiting... " << std::endl << endreq;
      exit(1);
    }

    return StatusCode::FAILURE;
  }
  std::string physFmt = "UNK";
  std::string fmtVersion;
  std::string dataIdent;

  // Get the information needed to find and interpret the bulk data:
  //   * physical storage type
  //   * version of the physical format
  //   * pathname or similar identifying information so the data can be found
  ret = m_meta->getReadInfo(ser, physFmt, fmtVersion, dataIdent);

  long storageType;
  StatusCode sc = decodeDescription(physFmt, storageType);

  // Depending on value of eDataFmt, figure out which private
  // conversion service to invoke.  Pass dataIdent and fmtVersion
  // as part of the opaque address.  
  // Create temporary address for the relevant type and classID 
  log << MSG::DEBUG << "Creating an address of type " 
      << (int)storageType << " for class " << classID << endreq;

  IOpaqueAddress* tmpAddress;
  const std::string par[3] = {dataIdent, fullpath, fmtVersion};
  const unsigned long ipar[1] = {ser};
  
  sc = addressCreator()->createAddress(storageType, classID, 
                                       par, ipar, tmpAddress);
  if ( !sc.isSuccess() ) {
    log << msgLevel 
	<< "Persistency service could not create a new address" << endreq;
    if (m_crash) {
      log << msgLevel << std::endl << "Exiting... " << std::endl << endreq;
      exit(1);
    }
     return sc;
  }  
  tmpAddress->addRef();

  // Set the transient store registry for the object associated to this address
  tmpAddress->setRegistry(entry);

  // Now create the object
  sc = m_detPersSvc->createObj(tmpAddress, refpObject);
  tmpAddress->release();
  if ( !sc.isSuccess() ) {
    log << msgLevel 
	<< "Persistency service could not create a new object" << endreq;
    if (m_crash) {
      log << msgLevel << std::endl << "Exiting... " << std::endl << endreq;
      exit(1);
    }
    return sc;
  }

  // Set validity of created object
  IValidity* pValidity = dynamic_cast<IValidity*>(refpObject);
  if (0 == pValidity) {
    log << MSG::WARNING
	<< "Created object does not implement IValidity: cannot set validity"
	<< endreq;
  } else {
    facilities::Timestamp* since;
    facilities::Timestamp* till;
    m_meta->getInterval(ser, since, till);
    log << MSG::DEBUG << "since, till " << since->getString() << " , " 
        << till->getString() << endreq;

    Gaudi::Time gSince = CalibData::CalibTime(*since).getGaudiTime();
    Gaudi::Time gTill = CalibData::CalibTime(*till).getGaudiTime();
    //log << MSG::DEBUG << "gSince " << gSince.year(false) << " " 
    //    << gSince.month(false) << " " << gSince.day(false) << " " 
    //    << gSince.hour(false) << " " << gSince.minute(false) << " " 
    //    <<gSince.second(false) << endreq;
    pValidity->setValidity(gSince, gTill);
  }

  log << MSG::DEBUG << "New object successfully created" << endreq;
  return StatusCode::SUCCESS;

}


/// Update a calib DataObject by calib type, flavor, time, and instrument
/// if necessary.
/// This method does not register DataObject in the transient data store,
/// The string storage type is discovered at runtime in the MySQL dbs.
/// Implementation:
/// - create a temporary address containing storage type and classID;
/// - dispatch to appropriate conversion service according to storage type;
/// - this will dispatch to appropriate converter according to CLID
///   (the CalibMySQLCnvSvc has no converters of its own).
StatusCode CalibMySQLCnvSvc::updateCalib( DataObject*        pObject,
                                          const std::string& fullpath,
                                          const Gaudi::Time&       time,
                                          const std::string& instr,
                                          const CLID&        classID,
                                          IRegistry*         entry )
{
  using CalibData::CalibBase;

  MsgStream log(msgSvc(), "CalibMySQLCnvSvc" );
  StatusCode status;

  // Look up calib object in the Metadata database
  std::string flavor = CalibData::CalibModelSvc::getFlavor(fullpath);
  std::string cType = CalibData::CalibModelSvc::getCalibType(fullpath);

  // ..and extra special munging for test
  if (std::string("Test") == cType.substr(0, 4)) {
    cType = std::string("Test_Gen");
  }

  if (0 == pObject) {
    log << MSG::ERROR << "There is no DataObject to update" << endreq;
    return StatusCode::FAILURE;
  }
  // Is object an instance of the specified class?
  if ( classID != pObject->clID() ) {
    log << MSG::ERROR << "Update requested for clID " << classID
	<< " while DataObject is of clID " 
	<< pObject->clID() << endreq;
    return StatusCode::FAILURE;
  }


  // check if already valid.  If so, we're done.
  //  Need access to IValidity interface
  CalibBase* pBase = dynamic_cast<CalibBase*>(pObject);
  if (pBase == 0) {
    log << MSG::ERROR 
        << "Object to be updated is not a calib object! " << endreq;
    return StatusCode::FAILURE;
  }

  if (pBase->isValid(m_detDataSvc->eventTime()) )
      return StatusCode::SUCCESS;

  // Following comes from createCalib.  Perhaps create and update
  // should be calling common utility since much of what they do is identical.
  unsigned int ser;
  calibUtil::Metadata::eRet ret;
  // Need to access db, so need connection
  m_metaConnected = true;

  if (m_useEventTime) {
    ret = m_meta->findBest(&ser, cType, CalibData::CalibTime(time),
                           m_qualityMask, instr, flavor);
  }
  else {
    ret = m_meta->findSoonAfter(&ser, cType, m_enterTimeStart, m_enterTimeEnd,
                                m_qualityMask, instr, flavor);
  }
  if (ret != calibUtil::Metadata::RETOk) {
    log << MSG::ERROR << "Could not access MySQL database" << endreq;
    return StatusCode::FAILURE;
  }
  else if (ser == 0) {  // no error, but no appropriate calib was found
    log << MSG::ERROR
        << "No appropriate calibration of (type, flavor) ("
        << cType << ", " << flavor << ")" << endreq;
    return StatusCode::FAILURE;
  }

  //  calibUtil::Metadata::eDataFmt physFmt = calibUtil::Metadata::FMTUnknown;
  std::string physFmt = "UNK";
  std::string fmtVersion;
  std::string dataIdent;

  // Get the information needed to find and interpret the bulk data:
  //   * physical storage type
  //   * version of the physical format
  //   * pathname or similar identifying information so the data can be found
  ret = m_meta->getReadInfo(ser, physFmt, fmtVersion, dataIdent);

  long storageType;
  status = decodeDescription(physFmt, storageType);

  // Depending on value of eDataFmt, figure out which private
  // conversion service to invoke.  Pass dataIdent and fmtVersion
  // as part of the opaque address.  
  // Create temporary address for the relevant type and classID 
  log << MSG::DEBUG << "Creating an address of type " 
      << (int)storageType << " for class " << classID << endreq;


  IOpaqueAddress* tmpAddress;
  const std::string par[3] = {dataIdent, fullpath, fmtVersion};
  const unsigned long ipar[1] = {ser};
  
  status = addressCreator()->createAddress(storageType, classID, 
                                           par, ipar, tmpAddress);
  if ( !status.isSuccess() ) {
    log << MSG::ERROR 
	<< "Persistency service could not create a new address" << endreq;
    return status;
  }  
  log << MSG::DEBUG << "Temporary address successfully created" << endreq;
  tmpAddress->addRef();

  // Set the transient store registry for the object associated to this address
  tmpAddress->setRegistry(entry);

  // Now update the object
  DataObject* pNewObject;
  status = m_detPersSvc->createObj(tmpAddress, pNewObject);
  tmpAddress->release();
  if ( !status.isSuccess() ) {
    log << MSG::ERROR 
	<< "Persistency service could not create object" << endreq;
    return status;
  }

  // Since DataObject::operator= operator is not virtual, dynamic cast first!
  // Overloaded virtual method Condition::update() must be properly defined!
  // The memory pointed to by the old pointer must contain the new object    
  // Note this dynamic cast only gets us part-way.  The object is
  // actually of some type derived from CalibBase.
  CalibBase* pNewBase = dynamic_cast<CalibBase*>(pNewObject);
  if  (0 == pNewBase) {
    log << MSG::ERROR
	<< "Cannot update objects other than Calib objects: " 
	<< "update() must be defined!"
	<< endreq;
    return StatusCode::FAILURE;
  }

  // Deep copy the new calib into the old DataObject.  To copy the *whole*
  // object, not just the CalibBase part, classes derived from CalibBase
  // must override update method.  
  // NOTE:  classes directly derived from CalibBase must call
  //        CalibBase::update in their own update method.
  pBase->update(*pNewBase, &log);  

  delete pNewBase;

  return StatusCode::SUCCESS;
}

StatusCode  CalibMySQLCnvSvc::decodeDescription(const std::string& description,
                                                long& type )
{
  MsgStream log(msgSvc(), "CalibMySQLCnvSvc");

  if (description == std::string("XML")) {
    type = XML_StorageType;
  }
  else if (description == std::string("ROOT")) {
    type = CALIBROOT_StorageType;
  }
  else {       // unsupported
    log << MSG::ERROR << "unsupported storage type " << description << endreq;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

/// Handle to the MySQL metadata database
calibUtil::Metadata* CalibMySQLCnvSvc::getMeta( ) {
  return m_meta;
}


StatusCode CalibMySQLCnvSvc::getValidInterval(unsigned int& serNo,
                                              Gaudi::Time** pvStart, Gaudi::Time** pvEnd) {
  using calibUtil::Metadata;
  using CalibData::CalibTime;

  if (*pvStart != 0) delete *pvStart;
  if (*pvEnd != 0) delete *pvEnd;
  facilities::Timestamp* since;
  facilities::Timestamp* till;
  Metadata::eRet ret = m_meta->getInterval(serNo, since, till);
  m_metaConnected = true;

  StatusCode status = StatusCode::FAILURE;

  if (ret == Metadata::RETOk) {
//    *pvStart = new CalibTime(*since);
//    *pvEnd = new CalibTime(*till);
      *pvStart = new Gaudi::Time(since->getClibTime(),since->getNano());
      *pvEnd = new Gaudi::Time(till->getClibTime(), till->getNano());
    status = StatusCode::SUCCESS;
  }

  delete since;
  delete till;
  return status;
}

/// At EndEvent check if MySQL connection is open periodically.  
/// If so, close it.
void CalibMySQLCnvSvc::handle(const Incident& inc) {
  static unsigned nEvent = 0;
  static unsigned freq = 5; 

  if (inc.type() != "EndEvent" ) return;

  nEvent++;
  if (m_metaConnected) {
    if  ((nEvent % freq) == 0) {
      m_meta->disconnectRead();
      m_metaConnected = false;
    }
  }

}

//void CalibMySQLCnvSvc::setCalibEnterTime(const ITime&  time, 
//                                         unsigned int interval)
// {
//  m_enterTime = time;
//  m_enterInterval = interval;
// }




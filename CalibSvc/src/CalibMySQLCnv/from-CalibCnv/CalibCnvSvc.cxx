// $Header$


// Do we really need all these?
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IConverter.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GenericAddress.h"

#include "CalibCnvSvc.h"
#include "facilities::Timestamp.h"
#include "calibUtil::Metadata.h"

// Instantiate our factory
static SvcFactory<CalibCnvSvc> calibcnvsvc_factory;
const ISvcFactory& CalibCnvSvcFactory = calibcnvsvc_factory;


// For the time being use GaudiKernel/ClassID.h predefined storage
// type MYSQL_StorageType, but would probably be better to have 
// user-defined type since our manner of using MySQL is probably 
// not what this type normally connotes.
CalibCnvSvc::CalibCnvSvc(const std::string& name, ISvcLocator* svc) :
  ConversionSvc(name, svc, MYSQL_StorageType), m_meta(0) {

  // Get our properties, if any.
  // Possibilities might include
  //    host for MySQL database
  //    table name


}

StatusCode CalibCnvSvc::initialize() {

  // Initialize ancestors
  StautsCode status = ConversionSvc::initialize();

  MsgStream log (msgSvc(), "CalibCnvSvc");

  if (!status.isSuccess() ) {
    return status;
  }

  if (!m_meta) m_meta = new Metadata();


  /* ++++ Not sure we want all this ++++ */

  // Locate the Calib Data Service.  Will use its service to
  // get timestamp current event.
  IDataProviderSvc* pProvider = 0;
  status = serviceLocator()->getService 
    ("CalibDataSvc",  IID_IDataProviderSvc, (IInterface*&)pProvider);
  if ( !status.isSuccess() ) {
    log << MSG::ERROR << "Could not locate CalibDataSvc" << endreq;
    return status;
  }

  // Set the CalibDataSvc as data provider service
  status = setDataProvider(pProvider);
  if ( !status.isSuccess() ) {
    log << MSG::ERROR << "Could not set data provider" << endreq;
    return status;
  }

  // Query the IDetDataSvc interface of the calib data service
  // Use it to get event time when creating addresses.
  // but since 
  status = pProvider->queryInterface(IID_IDetDataSvc, 
                                     (void**) &m_detDataSvc);
  if ( !status.isSuccess() ) {
    log << MSG::ERROR 
	<< "Cannot query IDetDataSvc interface of CalibDataSvc" 
	<< endreq;
    return status;
  } else {
    log << MSG::DEBUG 
	<< "Retrieved IDetDataSvc interface of CalibDataSvc" 
	<< endreq;
  }


  /* ++++ End of Not sure we want all this ++++ */


}

StatusCode CalibCnvSvc::finalize() {
  // disconnect from MySQL
  if (m_meta) {
    delete m_meta;
    m_meta = 0;
  }

  // delete parser instance

  // Shut down ancestors
  return ConversionSvc::finalize();
}


StatusCode CalibCnvSvc::queryInterface(const IID& riid, void** ppvInterface) {
  if (IID_ICalibCnvSvc.versionMatch(riid))  {
    *ppvInterface = static_cast<ICalibCnvSvc*> this;
  } else {
    // Interface is not directly availible: try out a base class
    return ConversionSvc::queryInterface(riid, ppvInterface);
  }
  addRef();
  return StatusCode::SUCCESS;
}


// This creates a very opaque address indeed: it saves enough information
// that the converter for the particular type of DataObject can 
// determine what file or files it will need to read:
// timestamp (kept as seconds since 1970)
// instrument type  (enumerated type or string??)
// calibration type (equivalent to class id, most likely;
//                   if so, don't need to keep as anonymous parameter)
// Now also need flavor -- use another std::string input for this?
//
//        For now
//            par[0] is instrument
//            par[1] is flavor
//
// Who calls this, anyway?
StatusCode createAddress(unsigned char svcType,
                         const CLID& clid,
                         const std::string* par,
                         const unsigned long* longPar,
                         IOpaqueAddress*& refpAddress) {
  if (MYSQL_StorageType != svcType) {
    log << MSG::ERROR
        << "Cannot create addresses of type " << (int)svcType
        << ", != "<< (int)MYSQL_StorageType  
        << endreq;
    return StatusCode::FAILURE;
  }

  // Get event time from somewhere (should be freom detDataSvc), convert 
  // to seconds.   For the time being (testing) just use "now"
  facilities::Timestamp now;
  int  timestamp = now.getClibTime();
  
  // Default to LAT.  Using enumerated type would be neater.
  // std::string instr = "LAT";

  refpAddress = new GenericAddress(MYSQL_StorageType,
                                   clid,
                                   par[0],  // string
                                   par[1],  // string
                                   timestamp);  // long int
  return StatusCode::SUCCESS;
}

// Implementation of ICalibCnvSvc interface
StatusCode CalibCnvSvc::findBest(std::string calibType,
                                 std::string instrument,
                                 facilities::Timestamp* eventTime,
                                 CalibKey *key) {
  using calibUtil::Metadata;

  eInstrument iInstr = translated(instrument);
  eRet status = m_meta->findBest(key, calibType, *eventTime, 
                                 LEVELProd || LEVELDev,
                                 iInstr);
}

StatusCode CalibCnvSvc::locatePDS(CalibKey key, std::string& dataFmt,
                                  std::string& fmtVersion,
                                  std::string& dataIdent) {

  eRet status;

  status = m_meta->getReadInfo(key, dataFmt, fmtVersion, dataIdent);

}

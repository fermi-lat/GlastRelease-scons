// $Header$

#include "GaudiKernel/IDetDataSvc.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/IConverter.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/GenericAddress.h"
#include "CalibRootCnvSvc.h"
#include "CalibData/CalibBase.h"
#include "cnv/RootBaseCnv.h"

// Make instances only via static factory class
static SvcFactory<CalibRootCnvSvc> calibRootCnvSvc_factory;
const ISvcFactory& CalibRootCnvSvcFactory = calibRootCnvSvc_factory;

CalibRootCnvSvc::CalibRootCnvSvc(const std::string& name, 
                               ISvcLocator* svc) :
  ConversionSvc(name, svc, CALIBROOT_StorageType),
  m_detPersSvc(0), m_detDataSvc(0)   {

  // Some day might have a property to declare having to do with path to
  // xml files.
}

StatusCode CalibRootCnvSvc::queryInterface(const InterfaceID& riid,
                                          void** ppvInterface) {
  /* Uncomment if choose to derive from abstract root conv. interface */
  if (IID_ICalibRootSvc.versionMatch(riid))  {
    *ppvInterface = (ICalibRootSvc*)this;
  }
  else {
    // Interface is not directly availible: try out a base class
    return ConversionSvc::queryInterface(riid, ppvInterface);
    /*  }  */
  addRef();
  }
  return StatusCode::SUCCESS;
}

StatusCode CalibRootCnvSvc::initialize() {
  StatusCode sc = ConversionSvc::initialize();

  MsgStream log(msgSvc(), "CalibRootCnvSvc");

  if (!sc.isSuccess()) return sc;

  // Locate the Calib Data Service.  Since it inherits from DataSvc
  // it has to implement IDataProviderSvc
  m_detDataSvc = 0;
  sc = serviceLocator()->getService 
    ("CalibDataSvc",  IID_IDataProviderSvc, (IInterface*&) m_detDataSvc);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR << "Could not locate CalibDataSvc" << endreq;
    return sc;
  }

  // Set the CalibDataSvc as data provider service
  sc = setDataProvider(m_detDataSvc);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR << "Could not set data provider" << endreq;
    return sc;
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
  sc = m_detPersSvc->queryInterface(IID_IAddressCreator, 
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
  log << MSG::DEBUG 
      << "Set it as the address creator of the CalibRootCnvSvc" << endreq;
  sc = setAddressCreator(iAddrCreator);
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR 	<< "Cannot set the address creator" << endreq;
    return sc;
  }

  // set properties if there are any??

  return sc;
}

StatusCode CalibRootCnvSvc::finalize() {
  // If anything was allocated, get rid of it.  So far, nothing was.

  return ConversionSvc::finalize();
}

StatusCode CalibRootCnvSvc::createAddress(long svc_type,
                                          const CLID& clid,
                                          const std::string* par, 
                                          const unsigned long* ip,
                                          IOpaqueAddress*& refpAddress) {

  MsgStream log( msgSvc(), name() );

  if ((unsigned) svc_type != CALIBROOT_StorageType) {
    log << MSG::ERROR << "bad storage type" << (int)svc_type << endreq;
    return StatusCode::FAILURE;
  }

  std::string dataIdent(par[0]); // file identifier for PDS version of data
  std::string fullpath(par[1]);  // path within TCDS for the object
  std::string fmtVersion(par[2]);
  int         serNo = ip[0];

  // for now have to ignore fmtVersion because of defective implementation
  // of GenericAddress. If we want it, should probably write new
  // opaque address implementation for this package to use.  All
  // dealings with (calibration) opaque addresses are confined to
  // the CalibSvc package.
  refpAddress = new GenericAddress(CALIBROOT_StorageType,
                                   clid,
                                   dataIdent,  
                                   fullpath,
                                   serNo);
                                   
  return StatusCode::SUCCESS; 

}

StatusCode CalibRootCnvSvc::writeToRoot(const std::string& outfile, 
                                        const std::string& tdsPath) {
  MsgStream log( msgSvc(), name() );

  // Find corresponding object
  DataObject* pObj;
  m_detDataSvc->findObject(tdsPath, pObj);
  if (!pObj) {
    log << "No object in TDS with path " << tdsPath << endreq;
    return StatusCode::FAILURE;
  }

  CalibData::CalibBase* pCalib = 
    dynamic_cast<CalibData::CalibBase*> (pObj);

  if (!pCalib) {
    log << "Object with path " << tdsPath << " not of proper type" << endreq;
    return StatusCode::FAILURE;
  }
  return writeToRoot(outfile, pCalib);
}
StatusCode CalibRootCnvSvc::writeToRoot(const std::string& outfile,
                                        CalibData::CalibBase* pCalib) {
  MsgStream log(msgSvc(), name() );

  // Find converter corresponding to this object
  IConverter* converter = ConversionSvc::converter(pCalib->clID());
  if (!converter) {
    log << "No converter found for object with CLID  " << pCalib->clID()
        << endreq;
    return StatusCode::FAILURE;
  }
  RootBaseCnv* pCnv = dynamic_cast<RootBaseCnv*>(converter);
  if (!pCnv) {
    log << "Converter for CLID " << pCalib->clID() <<  " not of proper type" 
        << endreq;
    return StatusCode::FAILURE;
  }
  // Call its createRoot method
  return pCnv->createRoot(outfile, pCalib);
}

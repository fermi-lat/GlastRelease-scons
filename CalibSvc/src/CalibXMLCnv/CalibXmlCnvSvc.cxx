// $Header$

// one or both of following includes may not be necessary..
// depends on whether we keep and use private members
// m_detPersSvc and m_detDataSvc

#include "GaudiKernel/IDetDataSvc.h"
#include "GaudiKernel/IConversionSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/GenericAddress.h"
#include "CalibXmlCnvSvc.h"
#include "xml/XmlParser.h"

// Make instances only via static factory class
static SvcFactory<CalibXmlCnvSvc> calibXmlCnvSvc_factory;
const ISvcFactory& CalibXmlCnvSvcFactory = calibXmlCnvSvc_factory;

CalibXmlCnvSvc::CalibXmlCnvSvc(const std::string& name, 
                               ISvcLocator* svc) :
  ConversionSvc(name, svc, XML_StorageType),
  m_detPersSvc(0), m_detDataSvc(0), m_parser(0) {

  // Some day might have a property to declare having to do with path to
  // xml files.
}

StatusCode CalibXmlCnvSvc::queryInterface(const IID& riid,
                                          void** ppvInterface) {
  if (IID_ICalibXmlSvc.versionMatch(riid))  {
    *ppvInterface = (ICalibXmlSvc*)this;
  }
  else {
    // Interface is not directly availible: try out a base class
    return ConversionSvc::queryInterface(riid, ppvInterface);
  }
  addRef();
  return StatusCode::SUCCESS;
}

StatusCode CalibXmlCnvSvc::initialize() {
  StatusCode status = ConversionSvc::initialize();

  MsgStream log(msgSvc(), "CalibXmlCnvSvc");

  if (!status.isSuccess()) return status;

  m_parser = new xml::XmlParser();

  if (m_parser == 0) {
    log << MSG::FATAL << "Unable to create XML parser instance \n"
        << endreq;
  }

  // set properties if there are any??

  return status;
}

StatusCode CalibXmlCnvSvc::finalize() {
  delete m_parser;
  m_parser = 0;
  return ConversionSvc::finalize();
}

StatusCode CalibXmlCnvSvc::createAddress(unsigned char svc_type,
                                         const CLID& clid,
                                         const std::string* par, 
                                         const unsigned long* ip,
                                         IOpaqueAddress*& refpAddress) {

  MsgStream log( msgSvc(), name() );

  if (svc_type != XML_StorageType) {
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
  refpAddress = new GenericAddress(XML_StorageType,
                                   clid,
                                   dataIdent,  
                                   fullpath,
                                   serNo);
                                   
  return StatusCode::SUCCESS; 

}




DOM_Document CalibXmlCnvSvc::parse(const char* filename) {
  return m_parser->parse(filename);
}

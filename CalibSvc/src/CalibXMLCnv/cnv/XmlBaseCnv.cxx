// $Header$

#include "XmlBaseCnv.h"

// #include <util/XMLUni.hpp>
// #include <util/XMLString.hpp>

#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/MsgStream.h"

#include "CalibSvc/ICalibXmlSvc.h"
#include "CalibSvc/ICalibMetaCnvSvc.h"

XmlBaseCnv::XmlBaseCnv( ISvcLocator* svc, const CLID& clid) :
  Converter (XML_StorageType, clid, svc),
  m_xmlSvc (0), m_metaSvc(0) {
}

StatusCode XmlBaseCnv::initialize() {
  StatusCode status = Converter::initialize();

  IDataProviderSvc* dp;

  // I guess the service names are assigned in jobOptions?

  serviceLocator()->getService ("CalibDataSvc",
                                IID_IDataProviderSvc,
                                (IInterface*&)dp);
  setDataProvider(dp);
  
  // Locate the Xml Conversion Service
  serviceLocator()->getService ("CalibXMLCnvSvc",
                                IID_ICalibXmlSvc,
                                (IInterface*&)m_xmlSvc);

  // Locate meta conversion service
  // Will anything need to be changed here to accommodate possibility
  // of two concrete implementations of ICalibMetaCnvSvc?  Would
  // have different storage types.  Could specify type desired
  // as job option
  serviceLocator()->getService("CalibMetaCnvSvc", 
                               IID_ICalibMetaCnvSvc,
                                (IInterface*&)m_metaSvc);
  return status;
}

StatusCode XmlBaseCnv::finalize() {
  return Converter::finalize();
}

// Create transient representation

StatusCode XmlBaseCnv::createObj(IOpaqueAddress* addr,
                                DataObject*&    refpObject) {

   // creates a msg stream for debug purposes
   MsgStream log( msgSvc(), "XmlBaseCnv" );
   
   if (0 == addr) {
     return StatusCode::FAILURE;
   }

  // first do the things we always need:
  //   First string parameter of opaque address is file ident
  //   Parse file into DOM representation
  const std::string* par = addr->par();
  DOM_Document doc = m_xmlSvc->parse(par[0].c_str());

  if (doc == DOM_Document() ) {
    log << MSG::FATAL 
        << "Unable to parse document " << par[0] << endreq;
    return StatusCode::FAILURE;
  }

  // Could conceivably write some code here to handle generic
  // parts of document.  Or, alternatively, add services to
  // CalibXmlCnvSvc for converters to invoke to do this.

  // Then do some fancy footwork in internalCreateObj to get the 
  // appropriate specific converter invoked to interpret the DOM 
  // correctly and make a new object of the correct kind.
  return internalCreateObj(doc, refpObject, addr);
}


/// In a backhanded way, invoke the right specific converter
/// for the type of the object to be created
StatusCode XmlBaseCnv::internalCreateObj(DOM_Document doc,
                                         DataObject*& refpObject,
                                         IOpaqueAddress* address) {
  // creates a msg stream for debug purposes
  MsgStream log( msgSvc(), "XmlBaseCnv" );
  
  // We're the default if we can't find anything better
  XmlBaseCnv* converter = this;      

  CLID classId = address->clID();

  IConverter* conv = this->conversionSvc()->converter(classId);

  if (0 == conv) {
    log << MSG::WARNING
        << "No proper converter found for classID " << classId
            << ", the default converter"
            << " will be used. " << endreq;
  } else {
    converter = dynamic_cast <XmlBaseCnv*> (conv);
    if (0 == converter) {
      log << MSG::ERROR
          << "The converter found for classID " << classId
              << " was not a descendent of XmlBaseCnv as it should be "
              << "( was of type " << typeid (*converter).name() << "). "
              << "The default converter will be used" << endreq;
      converter = this;
    }
  }
  
  // creates an object for the node found
  StatusCode sc = converter->i_createObj (doc, refpObject);
  if (sc.isFailure()) {
    return sc;
  }

  // ends up the object construction
  return converter->i_processObj(refpObject, address);
} 


// Default is to do nothing.  Derived classes may override.
StatusCode XmlBaseCnv::i_processObj(DataObject*, // pObject,
                                    IOpaqueAddress*) /* address) */ {
  return StatusCode::SUCCESS;
}
 

// Shouldn't ever really get here 
StatusCode XmlBaseCnv::i_createObj(DOM_Document, DataObject*&) {
  return StatusCode::FAILURE;
}






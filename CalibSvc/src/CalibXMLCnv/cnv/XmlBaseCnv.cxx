// $Header$

#include "XmlBaseCnv.h"

#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IAddressCreator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/MsgStream.h"

#include "CalibSvc/ICalibXmlSvc.h"
#include "CalibSvc/ICalibMetaCnvSvc.h"
#include "CalibData/CalibTime.h"
#include "CalibData/CalibBase.h"
#include "xml/Dom.h"

// A little ugly to include this here.  It's needed for 
// CAL-specific utilities findNextRange, findFirstRange
#include "idents/CalXtalId.h"

#include "facilities/Util.h"

#include <xercesc/dom/DOM_Document.hpp>
#include <xercesc/dom/DOM_NodeList.hpp>

// Local utilities to interpret attributes 
namespace {
  unsigned findRangeAtt(const DOM_Element& rangeElt) {
    using xml::Dom;
    using idents::CalXtalId;

    std::string att = Dom::getAttribute(rangeElt, "range");
    if (att.size() == 0) {  // range not applicable
      return 0;
    }
    if (att.compare(std::string("LEX8")) == 0) return CalXtalId::LEX8;
    if (att.compare(std::string("LEX1")) == 0) return CalXtalId::LEX1;
    if (att.compare(std::string("HEX8")) == 0) return CalXtalId::HEX8;
    if (att.compare(std::string("HEX1")) == 0) return CalXtalId::HEX1; 
    // anything else is illegal.  Should be caught by parser, but
    // maybe should also throw exception here.
  }

  unsigned findFace(const DOM_Element& faceElt) {
    using xml::Dom;
    using idents::CalXtalId;

    std::string att = Dom::getAttribute(faceElt, "end");
    if (att.compare(std::string("NEG")) == 0) return CalXtalId::NEG;
    if (att.compare(std::string("POS")) == 0) return CalXtalId::POS;
    return 0;     // in case "end" is not applicable, this is the answer
  }
}
    
XmlBaseCnv::~XmlBaseCnv() {}

// static CnvFactory<XmlBaseCnv> s_factory;
// const ICnvFactory& XmlBaseCnvFactory = s_factory;
XmlBaseCnv::XmlBaseCnv( ISvcLocator* svc, const CLID& clid) :
  Converter (XML_StorageType, clid, svc),
  m_xmlSvc (0), m_metaSvc(0), m_vstart(0), m_vend(0),
  m_nRow(10000), m_nCol(10000), m_nLayer(10000), m_nXtal(10000),
  m_nFace(10000), m_nRange(10000) {}

StatusCode XmlBaseCnv::initialize() {
  StatusCode status = Converter::initialize();

  IDataProviderSvc* dp;

  // I guess the service names are assigned in jobOptions?

  serviceLocator()->getService ("CalibDataSvc",
                                IID_IDataProviderSvc,
                                (IInterface*&)dp);
  setDataProvider(dp);
  
  // Locate the Xml Conversion Service
  serviceLocator()->getService ("CalibXmlCnvSvc",
                                IID_ICalibXmlSvc,
                                (IInterface*&)m_xmlSvc);

  // Locate meta conversion service
  // Will anything need to be changed here to accommodate possibility
  // of two concrete implementations of ICalibMetaCnvSvc?  Would
  // have different storage types.  Could specify type desired
  // as job option.  Ditto for name of class?
  serviceLocator()->getService("CalibMySQLCnvSvc", 
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

  std::string par0 = par[0];

  // Just in case there are environment variables in the file specification
  //  int nSub = 
  facilities::Util::expandEnvVar(&par0);

  //  DOM_Document doc = m_xmlSvc->parse(par[0].c_str());
  DOM_Document doc = m_xmlSvc->parse(par0.c_str());

  if (doc == DOM_Document() ) {
    log << MSG::FATAL 
        << "Unable to parse document " << par[0] << " aka " 
        << par0 << endreq;
    return StatusCode::FAILURE;
  }
  else {
    log << MSG::INFO
        << "successfully parsed document " << par[0] << " aka " 
        << par0 << endreq;
  }

  // Could conceivably write some code here to handle generic
  // parts of document.  Or, alternatively, add services to
  // CalibXmlCnvSvc for converters to invoke to do this.

  // Then do some fancy footwork in internalCreateObj to get the 
  // appropriate specific converter invoked to interpret the DOM 
  // correctly and make a new object of the correct kind.

  return internalCreateObj(doc.getDocumentElement(), refpObject, addr);
}


/** In a backhanded way, invoke the right specific converter
    for the type of the object to be created
    @param  elt      Document elt from XML document   (input)
    @param  refpObject 
*/
StatusCode XmlBaseCnv::internalCreateObj(const DOM_Element& docElt,
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

  unsigned int serNo = *(address->ipar());
  m_serNo = serNo;
  StatusCode sc = m_metaSvc->getValidInterval(serNo, 
                                              &m_vstart, 
                                              &m_vend );


  // creates an object for the node found
  if (sc.isSuccess()) sc = converter->i_createObj (docElt, refpObject);
  if (sc.isFailure()) {
    return sc;
  }

  // ends up the object construction
  sc = converter->i_processObj(refpObject, address);
  if (sc.isSuccess()) {
    log << MSG::DEBUG << "Successfully created calib. object " << endreq;
  }
  return sc;
} 

// Default is to do nothing.  Derived classes may override.
StatusCode XmlBaseCnv::i_processObj(DataObject*, // pObject,
                                    IOpaqueAddress* ) /* address */  {
  return StatusCode::SUCCESS;
}
 

// Shouldn't ever really get here 
StatusCode XmlBaseCnv::i_createObj(const DOM_Element&, DataObject*&) {
  return StatusCode::FAILURE;
}

/*
// Not sure yet whether this needs a real implementation or not
StatusCode XmlBaseCnv::updateObj(IOpaqueAddress* ,
                              DataObject*& ) {
  return StatusCode::FAILURE;
}

// Since we're not expecting to support writing back to persistent
// store, don't implement the converter *Rep functions.
*/

StatusCode XmlBaseCnv::readHeader(const DOM_Element&) {
  return StatusCode::SUCCESS;
}

const unsigned char XmlBaseCnv::storageType() {
  return XML_StorageType;
}

StatusCode XmlBaseCnv::readDimension(const DOM_Element& docElt, 
                                     unsigned& nRow, unsigned& nCol, 
                                     unsigned& nLayer,
                                     unsigned& nXtal, unsigned& nFace,
                                     unsigned& nRange) {
  using xml::Dom;

  DOM_Element dimElt = Dom::findFirstChildByName(docElt, "dimension");
  if (dimElt == DOM_Element()) return StatusCode::FAILURE;

  try {
    nRow = Dom::getIntAttribute(dimElt, "nRow");
    nCol = Dom::getIntAttribute(dimElt, "nCol");
    nLayer = Dom::getIntAttribute(dimElt, "nLayer");
    nXtal = Dom::getIntAttribute(dimElt, "nXtal");
    nFace = Dom::getIntAttribute(dimElt, "nFace");
    nRange = Dom::getIntAttribute(dimElt, "nRange");
  }
  catch (xml::DomException ex) {
    std::cerr << "From CalibSvc::XmlBaseCnv::readDimension" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }



  /*
  std::string att = Dom::getAttribute(dimElt, "nRow");
  nRow = (unsigned) atoi(att.c_str());
  att = Dom::getAttribute(dimElt, "nCol");
  nCol = (unsigned) atoi(att.c_str());
  att = Dom::getAttribute(dimElt, "nLayer");
  nLayer = (unsigned) atoi(att.c_str());
  att = Dom::getAttribute(dimElt, "nXtal");
  nXtal = (unsigned) atoi(att.c_str());
  att = Dom::getAttribute(dimElt, "nFace");
  nFace = (unsigned) atoi(att.c_str());
  att = Dom::getAttribute(dimElt, "nRange");
  nRange = (unsigned) atoi(att.c_str());
  */

  unsigned expected = nRow * nCol;
  // Make some consistency checks.  # tower elements should be nRow*nCol
  // # layer elements should be nRow*nCol*nLayer, etc.
  //  DOM_NodeList nlist = docElt.getElementsByTagName("tower");
  std::vector<DOM_Element> nlist;

  Dom::getDescendantsByTagName(docElt, "tower", nlist);
  
  if (nlist.size() != expected) {
    // put out a message and...
    return StatusCode::FAILURE;
  }
  expected *= nLayer;
  Dom::getDescendantsByTagName(docElt, "layer", nlist);
  if (nlist.size() != expected) {
    // put out a message and...
    return StatusCode::FAILURE;
  }

  expected *= nXtal;
  Dom::getDescendantsByTagName(docElt, "xtal", nlist);
  if (nlist.size() != expected) {
    // put out a message and...
    return StatusCode::FAILURE;
  }
  
  expected *= nFace;
  Dom::getDescendantsByTagName(docElt, "face", nlist);
  if (nlist.size() != expected) {
    // put out a message and...
    return StatusCode::FAILURE;
  }

  expected *= nRange;
  // Not as easy to check for the right number of range elements here
  // since they have different names, depending on what type of calibration
  // this is.  However they all have the same name, so just find first
  // child of a face element, then count all similarly-named elements

  DOM_Element child = Dom::getFirstChildElement(nlist[0]);
  std::string tagName = Dom::getTagName(child);
  Dom::getDescendantsByTagName(docElt, tagName, nlist);
  //  nlist = docElt.getElementsByTagName(child.getTagName());
  if (nlist.size() != expected) {
    // put out a message and...
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

/// Another utility for derived classes to use
void XmlBaseCnv::setBaseInfo(CalibData::CalibBase* pObj) {
  pObj->setValidity(*m_vstart, *m_vend);
  pObj->setSerNo(m_serNo);
}

DOM_Element XmlBaseCnv::findFirstRange(const DOM_Element& docElt) {
  using xml::Dom;
  using idents::CalXtalId;

  DOM_Element elt = Dom::findFirstChildByName(docElt, "tower");
  if (elt == DOM_Element()) return elt;

  //  std::string att = Dom::getAttribute(elt, "nRow");
  //  m_nRow = atoi(att.c_str());
  try {
    m_nRow = Dom::getIntAttribute(elt, "iRow");
    m_nCol = Dom::getIntAttribute(elt, "iCol");


    // All child elements of a tower are layer elements
    elt = Dom::getFirstChildElement(elt);
    //  att = Dom::getAttribute(elt, "iLayer");
    //  m_nLayer = atoi(att.c_str());
    m_nLayer = Dom::getIntAttribute(elt, "iLayer");

    // All child elements of layers are xtal elements
    elt = Dom::getFirstChildElement(elt);
    //    att = Dom::getAttribute(elt, "iXtal");
    //    m_nXtal = atoi(att.c_str());
    m_nXtal = Dom::getIntAttribute(elt, "iXtal");
  }
  catch (xml::DomException ex) {
    std::cerr << "From CalibSvc::XmlBaseCnv::findFirstRange" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }

  // All child elements of xtal are face elements
  elt = Dom::getFirstChildElement(elt); 
  m_nFace = findFace(elt);
  
  elt = Dom::getFirstChildElement(elt);
  m_nRange = findRangeAtt(elt);
  return elt;
}

/// Still another one to navigate XML file and find next set of range data
DOM_Element XmlBaseCnv::findNextRange(const DOM_Element& rangeElt) {
  using xml::Dom;

  DOM_Element elt = Dom::getSiblingElement(rangeElt);
  if (elt != DOM_Element()) {
    m_nRange = findRangeAtt(elt);
    return elt;
  }

  // Done with this xtal face; look for sibling
  DOM_Node node = rangeElt.getParentNode();
  elt = static_cast<DOM_Element &>(node);   // current xtal face
  elt = Dom::getSiblingElement(elt);          // next xtal face

  if (elt != DOM_Element()) {
    m_nFace = findFace(elt);

    elt = Dom::getFirstChildElement(elt);
    m_nRange = findRangeAtt(elt);
    return elt;
  }

  // Done with this xtal
  node = node.getParentNode();  // current xtal
  elt = static_cast<DOM_Element &>(node);
  elt = Dom::getSiblingElement(elt);         // next xtal

  if (elt != DOM_Element()) {
    try {
      m_nXtal = Dom::getIntAttribute(elt, "iXtal");
    }
    catch (xml::DomException ex) {
      std::cerr << "From CalibSvc::XmlBaseCnv::findNextRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    //    std::string att = Dom::getAttribute(elt, "iXtal");
    //    m_nXtal = atoi(att.c_str());

    // All child elements of xtal are face elements
    elt = Dom::getFirstChildElement(elt);
    m_nFace = findFace(elt);

    elt = Dom::getFirstChildElement(elt);
    m_nRange = findRangeAtt(elt);
    return elt;
  }

  // Done with this layer
  node = node.getParentNode();  // current layer
  elt = static_cast<DOM_Element &>(node);
  elt = Dom::getSiblingElement(elt);         // next layer

  if (elt != DOM_Element()) {
    //    std::string att = Dom::getAttribute(elt, "iLayer");
    //    m_nLayer = atoi(att.c_str());
    try {
      m_nLayer = Dom::getIntAttribute(elt, "iLayer");
    }
    catch (xml::DomException ex) {
      std::cerr << "From CalibSvc::XmlBaseCnv::findNextRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }
    

    // All child elements of layers are xtal elements
    elt = Dom::getFirstChildElement(elt);
    //    att = Dom::getAttribute(elt, "iXtal");
    //    m_nXtal = atoi(att.c_str());

    try {
      m_nXtal = Dom::getIntAttribute(elt, "iXtal");
    }
    catch (xml::DomException ex) {
      std::cerr << "From CalibSvc::XmlBaseCnv::findNextRange" << std::endl;
      std::cerr << ex.getMsg() << std::endl;
      throw ex;
    }

    // All child elements of xtal are face elements
    elt = Dom::getFirstChildElement(elt);
    m_nFace = findFace(elt);
    
    elt = Dom::getFirstChildElement(elt);
    m_nRange = findRangeAtt(elt);
    return elt;
  }

  // Done with this tower
  node = node.getParentNode();  // current tower
  elt = static_cast<DOM_Element &>(node);
  elt = Dom::getSiblingElement(elt);         // next tower

  if (elt == DOM_Element()) return elt;

  // otherwise we've got a new tower; go through the whole
  // rigamarole


  try {
    m_nRow = Dom::getIntAttribute(elt, "iRow");
    m_nCol = Dom::getIntAttribute(elt, "iCol");


    // All child elements of a tower are layer elements
    elt = Dom::getFirstChildElement(elt);
    //  att = Dom::getAttribute(elt, "iLayer");
    //  m_nLayer = atoi(att.c_str());
    m_nLayer = Dom::getIntAttribute(elt, "iLayer");

    // All child elements of layers are xtal elements
    elt = Dom::getFirstChildElement(elt);
    //    att = Dom::getAttribute(elt, "iXtal");
    //    m_nXtal = atoi(att.c_str());
    m_nXtal = Dom::getIntAttribute(elt, "iXtal");
  }
  catch (xml::DomException ex) {
    std::cerr << "From CalibSvc::XmlBaseCnv::findFirstRange" << std::endl;
    std::cerr << ex.getMsg() << std::endl;
    throw ex;
  }


  /*
  std::string att = Dom::getAttribute(elt, "iRow");
  m_nRow = atoi(att.c_str());

  att = Dom::getAttribute(elt, "iCol");
  m_nCol = atoi(att.c_str());

  // All child elements of a tower are layer elements
  elt = Dom::getFirstChildElement(elt);
  att = Dom::getAttribute(elt, "iLayer");
  m_nLayer = atoi(att.c_str());

  // All child elements of layers are xtal elements
  elt = Dom::getFirstChildElement(elt);
  att = Dom::getAttribute(elt, "iXtal");
  m_nXtal = atoi(att.c_str());
  */
  // All child elements of xtal are face elements
  elt = Dom::getFirstChildElement(elt);
  m_nFace = findFace(elt);

  elt = Dom::getFirstChildElement(elt);
  m_nRange = findRangeAtt(elt);

  return elt;
}

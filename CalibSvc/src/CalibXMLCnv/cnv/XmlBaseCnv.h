// $Header$
#ifndef CalibData_XmlBaseCnv_h
#define CalibData_XmlBaseCnv_h

/** @class XmlBaseCnv 

  Base class for calibration converters from XML files to TCDS.
  All such converters need to do certain things, which are
  handled here.

  @author J. Bogart
*/
#include "GaudiKernel/Converter.h"
#include <dom/DOM_Document.hpp>

class ISvcLocator;
class GenericAddress;
class ICalibXmlSvc;
class ICalibMetaSvc;

template <class TYPE> class CnvFactory;

class  XmlBaseCnv : public Converter {
public:
  virtual StatusCode initialize();

  virtual StatusCode finalize();

  /**
   Create the transient representation of an object, given an opaque
   address.  This and the following update method comprise the core 
   functionality of calibration converters.
  */
  virtual StatusCode createObj(IOpaqueAddress* addr,
                               DataObject*& refpObject);

  /**
     Update the transient representation of an object, given an opaque
     address.  This and the preceding create method comprise the core 
     functionality of calibration converters.
  */
  virtual StatusCode updateObj(IOpaqueAddress* addr,
                               DataObject*& refpObject);

  /** 
      Convert transient object to requested representation. This
      method is required since XmlBaseCnv inherits from Converter,
      but in practice it should not be called.
      Creation of persistent calibration data will 
      be accomplished by other means. 
  */
  virtual StatusCode createRep(DataObject*, // pObject,
                               IOpaqueAddress*&)  // refpAddress) 
  { return StatusCode::FAILURE;}

  /** 
      Updated converted representation of the transient object.  This
      method is required since XmlBaseCnv inherits from Converter,
      but in practice it should not be called.
  */
  virtual StatusCode updateRep(DataObject*,   // pObject,
                               IOpaqueAddress*&)  // refpAddress) 
  { return StatusCode::FAILURE; }

  ICalibXmlSvc* getCalibXmlSvc() {
    return m_xmlSvc;
  }

  static unsigned char storageType(){
    return XML_StorageType;
  }


protected:
  /**
     Constructor for this converter
     @param svc a ISvcLocator interface to find services
     @param clid the type of object the converter is able to convert
   */
  XmlBaseCnv(ISvcLocator* svc, const CLID& clid);

  virtual ~XmlBaseCnv() {};


  /** This creates the transient representation of an object from the
   *  DOM_Element representing it, then fills it and process it.
   *  This implementation actually only calls the i_* methods of the
   *  "right" converter to do the job; so the very first thing it
   *  does is get a pointer to the appropriate derived converter.
   *  Converters typically don't need to override this method
   *  but only to  override/implement some of the i_* methods.
   *  @param element the DOM_Element (actually a DOM_Document)
       to be used to builds the object
   *  @param refpObject the object to be built
   *  @param address the opaque address for this object
   *  @return status depending on the completion of the call
   */
  virtual StatusCode internalCreateObj (DOM_Document element,
                                        DataObject*& refpObject,
                                        IOpaqueAddress* address);
  
  /** This creates the transient representation of an object from the
   *  DOM_Element representing it. This actually does the "new" operation
   *  and deals with the attributes of the node. This base class implementation
   *  does nothing; it should not normally be called because it doesn't
   *  correspond to any TCDS class.  Instead, 
   *  i_createObj of some derived class will be called.
   *  @param element the DOM_Element (actually a DOM_Document) 
   *  to be used to builds the object
   *  @param refpObject the object to be built
   *  @return status depending on the completion of the call
   */
  virtual StatusCode i_createObj (DOM_Document element,
                                  DataObject*& refpObject);

  /// In case there is additional work to do on the created object
  virtual StatusCode i_processObj(DataObject* pObject,
                                  IOpaqueAddress* address);

  ICalibXmlSvc* m_xmlSvc;
  ICalibMetaSvc* m_metaSvc;

};


#endif

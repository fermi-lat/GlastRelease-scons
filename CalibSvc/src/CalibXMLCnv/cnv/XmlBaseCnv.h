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
#include "GaudiKernel/CnvFactory.h"
#include <xercesc/dom/DOM_Element.hpp>

class ISvcLocator;
class GenericAddress;
class ICalibXmlSvc;
class ICalibMetaCnvSvc;
class ITime;

namespace CalibData {
  class CalibTime;
  class CalibBase;
  class DacCol;
}

class  XmlBaseCnv : public Converter {

public:

  virtual ~XmlBaseCnv();

  virtual StatusCode initialize();

  virtual StatusCode finalize();

  /**
   Create the transient representation of an object, given an opaque
   address.  This and the following update method comprise the core 
   functionality of calibration converters.
  */
  virtual StatusCode createObj(IOpaqueAddress* addr,
                               DataObject*& refpObject);

  ICalibXmlSvc* getCalibXmlSvc() {
    return m_xmlSvc;
  }

  static const unsigned char storageType();

  /**
     Constructor for this converter
     @param svc a ISvcLocator interface to find services
     @param clid the type of object the converter is able to convert
   */
  //  XmlBaseCnv(ISvcLocator* svc, const CLID& clid = 0);
  XmlBaseCnv(ISvcLocator* svc, const CLID& clid);

protected:
  /** This creates the transient representation of an object from the
   *  DOM_Element representing it, then fills it and process it.
   *  This implementation actually only calls the i_* methods of the
   *  "right" converter to do the job; so the very first thing it
   *  does is get a pointer to the appropriate derived converter.
   *  Converters typically don't need to override this method
   *  but only to  override/implement some of the i_* methods.
   *  @param element the DOM_Element (typically the root element of the
   *   document) to be used to build the object
   *  @param refpObject the object to be built
   *  @param address the opaque address for this object
   *  @return status depending on the completion of the call
   */
  virtual StatusCode internalCreateObj (const DOM_Element& element,
                                        DataObject*& refpObject,
                                        IOpaqueAddress* address);
  
  /** This creates the transient representation of an object from the
   *  DOM_Element representing it. This actually does the "new" operation
   *  and deals with the attributes of the node. This base class implementation
   *  does nothing; it should not normally be called because it doesn't
   *  correspond to any TCDS class.  Instead, 
   *  i_createObj of some derived class will be called.
   *  @param element the DOM_Element (typically root element of document)
   *  to be used to builds the object
   *  @param refpObject the object to be built
   *  @return status depending on the completion of the call
   */
  virtual StatusCode i_createObj (const DOM_Element& element,
                                  DataObject*& refpObject);

  /// In case there is additional work to do on the created object
  virtual StatusCode i_processObj(DataObject* pObject,
                                  IOpaqueAddress* address);

  // Might want to verify that instrument, calType are correct,
  // for example.  If so, might as well provide the service in
  // the base converter.
  virtual StatusCode readHeader(const DOM_Element&);

  /// Convenience routine used by most CAL calibration types, which
  /// have a <dimension> element describing how the remainder of the
  /// data is laid out.
  StatusCode readDimension(const DOM_Element& docElt, 
                           unsigned& nRow, unsigned& nCol, 
                           unsigned& nLayer,
                           unsigned& nXtal, unsigned& nFace,
                           unsigned& nRange,
                           unsigned* nDacCol=0);

  /// Another one to find first range element
  DOM_Element findFirstRange(const DOM_Element& docElt);

  /// Still another one to navigate XML file and find next set of range data
  DOM_Element findNextRange(const DOM_Element& rangeElt);

  /// Another one to find first dac collection element
  DOM_Element findFirstDacCol(const DOM_Element& docElt);

  /// Still another one to navigate XML file and find next dac collection
  DOM_Element findNextDacCol(const DOM_Element& rangeElt);

  CalibData::DacCol* processDacCol(DOM_Element dacColElt, unsigned* range);

  /// Another convenience for derived classes: sets information belonging
  /// to the calibration base class, namely validity interval and serial
  /// number.
  void setBaseInfo(CalibData::CalibBase* pObj);

  ICalibXmlSvc* m_xmlSvc;
  ICalibMetaCnvSvc* m_metaSvc;

  int m_serNo;
  ITime*  m_vstart;
  ITime*  m_vend;

  /// A place to keep track of where we are if we're handling CAL data
  unsigned m_nRow;
  unsigned m_nCol;
  unsigned m_nLayer;
  unsigned m_nXtal;
  unsigned m_nFace;
  unsigned m_nRange;

};

#endif

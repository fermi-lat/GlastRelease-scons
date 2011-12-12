// $Header$
#ifndef CalibData_XmlBaseCnv_h
#define CalibData_XmlBaseCnv_h

/** @class XmlBaseCnv 

  Base class for calibration converters from XML files to TCDS.
  All such converters need to do certain things, which are
  handled here.

  @author J. Bogart
*/
#include <string>
#include <vector>
#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/Time.h"
#include <xercesc/dom/DOMElement.hpp>

class ISvcLocator;
class GenericAddress;
class ICalibXmlSvc;
class ICalibMetaCnvSvc;
//class ITime;

namespace CalibData {
  class CalibTime;
  class CalibBase;
  class DacCol; // for now used only by calorimeter intNonLin calibration
  class Xpos;   // only of interest for some calorimeter calibrations
  class ValSig;
}

using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;


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
   *  DOMElement representing it, then fills it and process it.
   *  This implementation actually only calls the i_* methods of the
   *  "right" converter to do the job; so the very first thing it
   *  does is get a pointer to the appropriate derived converter.
   *  Converters typically don't need to override this method
   *  but only to  override/implement some of the i_* methods.
   *  @param element the DOMElement (typically the root element of the
   *   document) to be used to build the object
   *  @param refpObject the object to be built
   *  @param address the opaque address for this object
   *  @return status depending on the completion of the call
   */
  virtual StatusCode internalCreateObj (const DOMElement* element,
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
  virtual StatusCode i_createObj (const DOMElement* element,
                                  DataObject*& refpObject);

  /// In case there is additional work to do on the created object
  virtual StatusCode i_processObj(DataObject* pObject,
                                  IOpaqueAddress* address);

  // Might want to verify that instrument, calType are correct,
  // for example.  If so, might as well provide the service in
  // the base converter.
  virtual StatusCode readHeader(const DOMElement*);

  /// Retrieve the class type of the data store the converter uses.
  virtual long repSvcType() const {return Converter::i_repSvcType();}


  /// Find first range element.  Derived classes which need it
  /// must define their own implementation.
  DOMElement* findFirstRange(const DOMElement* /* docElt */) {
    return 0;}


  /// Still another one to navigate XML file and find next set of range data
  DOMElement* findNextRange(const DOMElement* /* rangeElt */) {
    return 0;}

  /// Another one to find first dac collection element
  DOMElement* findFirstDacCol(const DOMElement* docElt);

  /// Still another one to navigate XML file and find next dac collection
  DOMElement* findNextDacCol(const DOMElement* rangeElt);

  CalibData::DacCol* processDacCol(DOMElement* dacColElt, unsigned* range);

  DOMElement* findXpos(const DOMElement* docElt);

  CalibData::Xpos* processXpos(DOMElement* xposElt);

  /// Read in what will become a CalibData::ValSig
  CalibData::ValSig* processValSig(DOMElement* elt, 
                                   std::string valName, std::string sigName);

  /// Read in what will become a vector of CalibData::ValSig
  std::vector<CalibData::ValSig>*  processValSigs(DOMElement* elt, 
                                                  std::string valName, 
                                                  std::string sigName);

  /// Another convenience for derived classes: sets information belonging
  /// to the calibration base class, namely validity interval and serial
  /// number.
  void setBaseInfo(CalibData::CalibBase* pObj);

  ICalibXmlSvc* m_xmlSvc;
  ICalibMetaCnvSvc* m_metaSvc;

  int m_serNo;
  Gaudi::Time*  m_vstart;
  Gaudi::Time*  m_vend;


};

#endif

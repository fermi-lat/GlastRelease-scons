//$Header$
#ifndef CalibXmlCnvSvc_h
#define CalibXmlCnvSvc_h  1

#include <string>

#include "CalibSvc/ICalibMetaCnvSvc.h"   // not sure if this is needed
//   #include "CalibSvc/IInstrumentName.h"
#include "GaudiKernel/ConversionSvc.h"
#include "CalibSvc/ICalibXmlSvc.h"

#include <dom/DOM_Document.hpp>

/// Forward and external declarations
template <class TYPE> class SvcFactory;

class IDetDataSvc;
class IOpaqueAddress;
namespace xml {
  class XmlParser;
}

///---------------------------------------------------------------------------
/** @class CalibXmlCnvSvc

    A conversion service for GLAST calibration bulk data in Xml format.

    @author J. Bogart
    @date November 2002
*///--------------------------------------------------------------------------

class CalibXmlCnvSvc : public ConversionSvc, virtual public ICalibXmlSvc 
{
  /// Only factories can access protected constructors
  friend class SvcFactory<CalibXmlCnvSvc>;

 protected:

  CalibXmlCnvSvc(const std::string& name, ISvcLocator* svc );
  virtual ~CalibXmlCnvSvc() {}

 public:
  
  // Reimplemented from IInterface

  virtual StatusCode queryInterface( const IID& riid, 
				     void** ppvInterface);  

 public:

  // Overloaded from ConversionSvc

  virtual StatusCode initialize();
  virtual StatusCode finalize();

  /**
   * Create an XML address using explicit arguments to identify a single object
   * @param svc_type the service type
   * @param CLID the CLID of the XML Element for which an address is created
   * @param par an array of three strings containing the format version,
   *        calibration type name and the flavor, in this order
   * @param ip has a single element, the serial number of the MySQL row
   *        which corresponds to this element
   * @param refpAddress the new address created
   * @return a StatusCode giving the status of the address creation
   */
  virtual StatusCode createAddress(unsigned char svc_type,
                                   const CLID& clid,
                                   const std::string* par, 
                                   const unsigned long* ip,
                                   IOpaqueAddress*& refpAddress);

  // from ICalibXmlSvc interface
  virtual DOM_Document parse(const char* filename);

  
  /*  
         probably don't need to override these
  /// Create a transient representation from another rep of this object.
  virtual StatusCode createObj     ( IOpaqueAddress* pAddress, 
				     DataObject*&    refpObject );
  
  /// Resolve the references of the created transient object.
  virtual StatusCode fillObjRefs   ( IOpaqueAddress* pAddress, 
				     DataObject* pObject );
  
  /// Update a transient representation from another rep of this object.
  virtual StatusCode updateObj     ( IOpaqueAddress* pAddress, 
				     DataObject* pObject );

  /// Update the references of an updated transient object.
  virtual StatusCode updateObjRefs ( IOpaqueAddress* pAddress, 
				     DataObject* pObject );

  /// Convert a transient object to a requested representation.
  virtual StatusCode createRep     ( DataObject* pObject, 
				     IOpaqueAddress*& refpAddress );

  /// Resolve the references of a converted object. 
  virtual StatusCode fillRepRefs   ( IOpaqueAddress* pAddress,
				     DataObject* pObject );

  /// Update a converted representation of a transient object.
  virtual StatusCode updateRep     ( IOpaqueAddress* pAddress, 
				     DataObject* pObject );

  /// Update the references of an already converted object.
  virtual StatusCode updateRepRefs ( IOpaqueAddress* pAddress, 
				     DataObject* pObject );
  */


  
 public:

  // Probably don't need this.  meta conversion service should
  // already have fetched and stored everything we need.
  //  calibUtil::Metadata* getMeta();

 private:

  //Probably don't need this
  // /// Handle for metadata access
  // calibUtil::Metadata*    m_meta;

  // Not sure we need this
  /// Handle to the IConversionSvc interface of the DetectorPersistencySvc
  IConversionSvc*      m_detPersSvc;

  // Not sure we need this
  /// Handle to the IDetDataSvc interface of the CalibDataSvc
  IDetDataSvc*         m_detDataSvc;

  xml::XmlParser*       m_parser;
};
#endif   









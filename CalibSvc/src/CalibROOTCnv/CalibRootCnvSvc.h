//$Header$
#ifndef CalibRootCnvSvc_h
#define CalibRootCnvSvc_h  1

#include <string>

#include "GaudiKernel/ConversionSvc.h"


//  Hi Heather.  For the XML conversion service I defined an
//  extra interface to do generic XML things that converters
//  might need.  For now there really is just one thing there:
//  a function which reads in the XML file, so that the converters
//  deal with the DOM representation (in-memory, but otherwise 
//  isomorphic to the XML description) rather than with the physical
//  file.  Putting that kind of functionality in the conversion service
//  is probably a good idea, but it doesn't have to be defined in
//  an abstract interface.  In retrospect, I'd be tempted to just 
//  let the converters have direct access to the conversion service
//  implementation and not bother with ICalibXmlSvc.
//  If there are similar things that the ROOT conversion service 
//  might do, they belong inside CalibRootCnvSvc as public methods
//  (and also in an abstract interface if you decide you want one).

// #include "CalibSvc/ICalibXmlSvc.h"

/// Forward and external declarations
template <class TYPE> class SvcFactory;

class IDetDataSvc;
class IOpaqueAddress;


///---------------------------------------------------------------------------
/** @class CalibRootCnvSvc

    A conversion service for GLAST calibration bulk data in ROOT format.

    @author H Kelly, J. Bogart
    @date February 2003
*///--------------------------------------------------------------------------

class CalibRootCnvSvc : public ConversionSvc
{
  /// Only factories can access protected constructors
  friend class SvcFactory<CalibRootCnvSvc>;

 protected:

  CalibRootCnvSvc(const std::string& name, ISvcLocator* svc );
  virtual ~CalibRootCnvSvc() {}

 public:
  
  // Reimplemented from IInterface

  virtual StatusCode queryInterface( const IID& riid, 
				     void** ppvInterface);  

 public:

  // Overloaded from ConversionSvc

  virtual StatusCode initialize();
  virtual StatusCode finalize();

  /**
   * Create a ROOT address using explicit arguments to identify a single object
   * @param svc_type the service type
   * @param CLID the CLID of the ROOT Element for which an address is created
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


  /*  
      There are a pile of functions implemented in the Gaudi
      base class ConversionSvc which we can just let be, such
      as createObj, fillObjRefs,...

      The base implementation looks up the appropriate converter
      and invokes it, usually just what we want to do.
  */
 private:

  // With current functionality, there really is no need to keep
  // this as a member.  It's only used during initialize()
  /// Handle to the IConversionSvc interface of the DetectorPersistencySvc
  IConversionSvc*      m_detPersSvc;

  // Ditto (don't really need to keep this either)
  /// Handle to the IDetDataSvc interface of the CalibDataSvc
  IDetDataSvc*         m_detDataSvc;

};
#endif   









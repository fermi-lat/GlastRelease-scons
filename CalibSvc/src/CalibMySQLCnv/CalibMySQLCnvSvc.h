//$Header$
#ifndef CalibMySQLCnvSvc_h
#define CalibMySQLCnvSvc_h  1

/// Include files
#include "CalibSvc/ICalibMetaCnvSvc.h"
#include "GaudiKernel/ConversionSvc.h"

/// Forward and external declarations
// class ConditionsDBGate;
template <class TYPE> class SvcFactory;
class IDetDataSvc;
class IOpaqueAddress;

///---------------------------------------------------------------------------
/** @class CalibMySQLCnvSvc

    A conversion service for GLAST calibration metadata database persistency.
    Allows to create and update condition data objects (i.e. DataObjects
    implementing IValidity).

    Adapted from LHCb class DetCond/ConditionsDBCnvSvc by Andrea Valassi
    @author J. Bogart
    @date November 2002
*///--------------------------------------------------------------------------

class CalibMySQLCnvSvc : public ConversionSvc, 
                         virtual public ICalibMetaCnvSvc
{
  /// Only factories can access protected constructors
  friend class SvcFactory<CalibMySQLCnvSvc>;

 protected:

  CalibMySQLCnvSvc(const std::string& name, ISvcLocator* svc );
  virtual ~CalibMySQLCnvSvc();

 public:
  
  // Reimplemented from IInterface

  virtual StatusCode queryInterface( const IID& riid, 
				     void** ppvInterface );  

 public:

  // Overloaded from ConversionSvc

  virtual StatusCode initialize();
    virtual StatusCode finalize();
  
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

  /// Create an address using explicit arguments to identify a single object.
  virtual StatusCode createAddress ( unsigned char svc_type,
				     const CLID& clid,
				     const std::string* par, 
				     const unsigned long* ip,
				     IOpaqueAddress*& refpAddress );
  
 public:
  // Implementation of IConditionsDBCnvSvc.
  // Create/update condition DataObject not necessarily registered in the TDS.
  
  /// Create a condition DataObject by calib type name, flavor and time.
  /// This method does not register DataObject in the transient data store,
  /// but may register TDS addresses for its children if needed (e.g. Catalog).
  /// The string storage type is discovered at runtime in the Metadata dbs.
  /// The entry name identifies a condition amongst the many in the string.
  StatusCode createCalibMetadata(DataObject*& refpObject,
                                 const std::string& calibType,
                                 const std::string& flavor,
                                 const ITime&       time,
                                 const std::string& instrument,
                                 const CLID&        classID,
                                 IRegistry*         entry=0) = 0;
  
  /// Update a condition DataObject by 
  /// This method does not register DataObject in the transient data store,
  /// but may register TDS addresses for its children if needed (e.g. Catalog).
  /// The string storage type is discovered at runtime in the CondDB.
  /// The entry name identifies a condition amongst the many in the string.
  StatusCode updateCalibMetadata (DataObject*          pObject,
                                 const std::string& calibType,
                                 const std::string& flavor,
                                 const ITime&       time,
                                 const CLID&        classID,
                                 IRegistry*         entry=0) = 0;
  
  /// Decode, decode the string storage type from/to the folder descrip string
  //  StatusCode decodeDescription     ( const std::string&   description,
  //				     unsigned char&       type);
  //  StatusCode encodeDescription     ( const unsigned char& type,
  //			     std::string&         description);
  // Above just used ascii rep. of storage type in dbs; our dbs has strings
  // so need a way to associate them with the storage types.  
  /// Handle to metadata connection
  calibUtil::Metadata* meta();

 private:

  /// Handle for metadata access
  calibUtil::Metadata*    m_meta;

  /// Handle to the IConversionSvc interface of the DetectorPersistencySvc
  IConversionSvc*      m_detPersSvc;

  /// Handle to the IDetDataSvc interface of the CalibDataSvc
  IDetDataSvc*         m_detDataSvc;

  // Maybe also handle to IInstrumentName interface? Or do we really
  // need all these different handles to the same class??
  IInstrumentNameSvc*  m_instrSvc;
};
#endif   


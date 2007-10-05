//$Header$
#ifndef MootCnvSvc_h
#define MootCnvSvc_h  1

/// Include files
#include <vector>
#include <string>
#include <cstring>

#ifdef _WIN32
#  include <hash_map>
using  stdext::hash_map;
// using  stdext::hash;
#else
#  include <ext/hash_map>
using __gnu_cxx::hash_map;
// using __gnu_cxx::hash;
#endif

#include "CalibSvc/IMootSvc.h"
#include "CalibSvc/IInstrumentName.h"
#include "GaudiKernel/ConversionSvc.h"
#include "GaudiKernel/IDataManagerSvc.h"

/// Forward and external declarations
// class ConditionsDBGate;
template <class TYPE> class SvcFactory;
class IDetDataSvc;
class IOpaqueAddress;
class MsgStream;



namespace MOOT {
  class MoodConnection;
  class MootQuery;
}


// Just make up something unlikely to be used elsewhere
const long MOOT_StorageType = POOL_StorageType - 1;


///---------------------------------------------------------------------------
/** @class MootSvc

    A conversion service for GLAST calibration metadata database persistency.
    Allows to create and update condition data objects (i.e. DataObjects
    implementing IValidity).

    Adapted from LHCb class DetCond/ConditionsDBCnvSvc by Andrea Valassi
    @author J. Bogart
    @date November 2002
*///--------------------------------------------------------------------------

class MootSvc : public ConversionSvc, virtual public IMootSvc
{
  /// Only factories can access protected constructors
  friend class SvcFactory<MootSvc>;

#ifdef _WIN32
typedef hash_map<const char*, int> HashMap; 
//                 hash<const char*, eqstr> > HashMap;
#else
typedef hash_map<const char*, int> HashMap;
                 //                 hash<const char*, eqstr> > HashMap;
#endif

  typedef std::pair<const char*, int>  hashpair;
 public:

  // Stuff for application clients  (implementing IMootSvc interface)
  virtual std::string getLatcSourcePath()  {return m_latcPath;}

  // Reimplemented from IInterface

  virtual StatusCode queryInterface( const InterfaceID& riid, 
				     void** ppvInterface );  

  virtual StatusCode initialize();
  virtual StatusCode finalize();

  // stuff for other services within the package, like Moot converters
  unsigned getLatcParmMaxCnt();


 protected:

  MootSvc(const std::string& name, ISvcLocator* svc );
  virtual ~MootSvc();

  // Overloaded from ConversionSvc
  // Suspect we don't need any of them.
  /*  
  /// Create a transient representation from another rep of this object.
  virtual StatusCode createObj     ( IOpaqueAddress* pAddress, 
				     DataObject*&    refpObject );
  
  // Maybe don't need to override
  
  /// Update a transient representation from another rep of this object.
  virtual StatusCode updateObj     ( IOpaqueAddress* pAddress, 
				     DataObject* pObject );

  /// Create an address using explicit arguments to identify a single object.
  virtual StatusCode createAddress ( long svc_type,
				     const CLID& clid,
				     const std::string* par, 
				     const unsigned long* ip,
				     IOpaqueAddress*& refpAddress );
  */  
 public:
  // Implementation of IMootSvc.
  virtual int latcParmIx(const std::string& parmClass);
 

  // For internal CalibSvc use, but not necessarily MootSvc class
  virtual MOOT::MootQuery* getConnection(bool verbose=true);

  // To be called only by CalibDataSvc, so it's public but is not in
  // the public include directory and is not part of the IMootSvc interface.
  // For now we handle only parameter classes belonging to a precinct
  // other than 'generic'
  StatusCode makeMootNodes(const std::string& parent);

  

 private:

  StatusCode getPrecincts();  // fills m_prNames

  /// Handles for metadata access
  MOOT::MootQuery*    m_q;
  MOOT::MoodConnection* m_c;

  /// MOOT archive path.  
  std::string   m_archive;  

  // Precinct names.  Index comes from public enum
  std::vector<std::string> m_prNames;
  /*
    probably should have a job option.  Special value of * would mean
    'use existing value of env var MOOT_ARCHIVE'   Otherwise 
    do setenv of MOOT_ARCHIVE, using supplied value
  */


  /// Handle to the IConversionSvc interface of the DetectorPersistencySvc
  IConversionSvc*      m_detPersSvc;

  /// Handle to the IDataProviderSvc interface of the CalibDataSvc
  IDataProviderSvc*    m_provider;

  /// Handle to the IDataProviderSvc interface of the CalibDataSvc
  IDataManagerSvc*    m_dataManager;

  //  IAddressCreator*     m_addrCreator;

  IInstrumentName*     m_instrSvc;   // might not need this one

  MsgStream*           m_log;

  unsigned             m_master; // latc master from most recent event
  std::string          m_latcPath;
  //  HashMap*             m_latcParmMap;
  hash_map<const char*, int>* m_latcParmMap;
};
#endif   


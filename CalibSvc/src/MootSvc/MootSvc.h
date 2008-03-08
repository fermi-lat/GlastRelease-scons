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

class IDataProviderSvc;

#include "GaudiKernel/Service.h"
#include "CalibSvc/IMootSvc.h"
#include "CalibSvc/IInstrumentName.h"


/// Forward and external declarations
// class ConditionsDBGate;
template <class TYPE> class SvcFactory;
class ISvcLocator;
class MsgStream;

namespace MOOT {
  class MoodConnection;
  class MootQuery;
}

namespace CalibData {
  class MootParmCol;
}

/** @class MootSvc
    Implements IMootSvc interface, providing access to Moot data.
    See also CalibData/Moot for definition of data structures
    
    @author J. Bogart
*/
class MootSvc :  public Service, virtual public IMootSvc
{
  /// Only factories can access protected constructors
  friend class SvcFactory<MootSvc>;

#ifdef _WIN32
typedef hash_map<const char*, int> HashMap; 

#else
typedef hash_map<const char*, int> HashMap;

#endif

  typedef std::pair<const char*, int>  hashpair;
 public:

  // Reimplemented from IMootSv
  virtual const CalibData::MootParmCol* getMootParmCol(unsigned& hw);

  /// Return last LATC master key seen in data
  virtual unsigned getHardwareKey() ;


  /// Return index in MootParmCol of specified class
  virtual int latcParmIx(const std::string& parmClass) const;

  virtual MOOT::MootQuery* getConnection() const {return m_q;}

  // Reimplemented from IInterface

  virtual StatusCode queryInterface( const InterfaceID& riid, 
				     void** ppvInterface );  


  virtual StatusCode initialize();
  virtual StatusCode finalize();

  // stuff for other services within the package,
  unsigned getLatcParmMaxCnt();


 protected:

  MootSvc(const std::string& name, ISvcLocator* svc );
  virtual ~MootSvc();

 public:

  // For internal CalibSvc use, but not necessarily MootSvc class
  // Not part of IMootSvc interface
  virtual MOOT::MootQuery* makeConnection(bool verbose=true);
  typedef std::pair<std::string, std::string> ParmPrecinct;
  const std::vector<ParmPrecinct>* getParmPrecinct() const 
  {return  &m_parmPrecincts;}

 private:

  StatusCode getPrecincts();  // fills m_prNames, hash map
  StatusCode updateMootParmCol();
  StatusCode updateFswKeys();

  /// Handles for metadata access
  MOOT::MootQuery*    m_q;
  MOOT::MoodConnection* m_c;

  /// MOOT archive path.  
  std::string   m_archive;  

  // Precinct names.  Index comes from public enum
  std::vector<std::string> m_prNames;
  hash_map<const char*, int>* m_parmMap;

  std::vector<ParmPrecinct> m_parmPrecincts;

  /*
    probably should have a job option.  Special value of * would mean
    'use existing value of env var MOOT_ARCHIVE'   Otherwise 
    do setenv of MOOT_ARCHIVE, using supplied value
  */



  IInstrumentName*     m_instrSvc;   // might not need this one
  IDataProviderSvc* m_eventSvc;

  MsgStream*           m_log;

  unsigned             m_hw; // latc master from most recent event
  //  std::string          m_latcPath;

  bool  m_useEventKeys;                       // job options

  CalibData::MootParmCol*  m_mootParmCol;
};
#endif   


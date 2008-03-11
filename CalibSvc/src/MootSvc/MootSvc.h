//$Header$
#ifndef MootCnvSvc_h
#define MootCnvSvc_h  1

/// Include files
#include <vector>
#include <string>
#include <cstring>

class IDataProviderSvc;

#include "GaudiKernel/Service.h"
#include "CalibSvc/IMootSvc.h"
#include "CalibSvc/IInstrumentName.h"


/// Forward and external declarations
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

 public:

  // Reimplemented from IMootSv

  /// Return absolute path for parameter source file of specified class.
  /// If none return empty string.
  std::string getMootParmPath(const std::string& cl, unsigned& hw);

  /// Return MootParm structure for parameter source file of specified class.
  /// If none return blank structure.
  virtual const CalibData::MootParm* getMootParm(const std::string& cl, 
                                           unsigned& hw);

  // Return pointer to Moot parameter collection.  Also set output
  // arg. hw to current hw key
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

 protected:

  MootSvc(const std::string& name, ISvcLocator* svc );
  virtual ~MootSvc();

 public:

  // For internal CalibSvc use, but not necessarily MootSvc class
  // Not part of IMootSvc interface
  virtual MOOT::MootQuery* makeConnection(bool verbose=true);
  typedef std::pair<std::string, std::string> ParmPrecinct;


 private:

  StatusCode getPrecincts();  // fills m_parmPrecincts

  /// Given param class name, find associated precinct
  std::string findPrecinct(const std::string& pclass);

  StatusCode updateMootParmCol();
  StatusCode updateFswKeys();

  /// Handles for metadata access
  MOOT::MootQuery*    m_q;
  MOOT::MoodConnection* m_c;

  /// MOOT archive path.  
  std::string   m_archive;  

  std::vector<ParmPrecinct> m_parmPrecincts;

  IInstrumentName*     m_instrSvc;   // might not need this one
  IDataProviderSvc* m_eventSvc;

  MsgStream*           m_log;

  unsigned             m_hw; // latc master from most recent event

  bool  m_useEventKeys;                       // job options
  bool  m_verbose;                            // controls MoodConnection

  CalibData::MootParmCol*  m_mootParmCol;
};
#endif   


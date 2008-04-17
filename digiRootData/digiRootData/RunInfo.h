#ifndef ROOT_RUNINFO_H
#define ROOT_RUNINFO_H 1

#include "Riostream.h"
#include "TObject.h"
#include "enums/Lsf.h"

/** @class RunInfo
* @brief Encapsulate information about the run that an event was taken in
*
*        The RunInfo consists of:
*            platform    -> the type of platfrom that made the data ( One of enums::Lsf::Platform )
*            dataOrigin  -> to distinguish ground, flight and MC ( One of enums::Lsf::DataOrigin )
*            id          -> run number assigned on the ground. 
*            startTime   -> number of seconds since start of GLAST epoch at which the run started
*            dataTransferId -> number of seconds since start of GLAST epoch
*            at which the run started
*
* $Header$
*/

class RunInfo : public TObject {
    
public:
  
  /// Default c'tor.  Assigns sentinel values to all fields
  RunInfo( )
    : m_platform(enums::Lsf::NoPlatform), m_origin(enums::Lsf::NoOrigin), 
      m_id(LSF_INVALID_UINT), m_startTime(LSF_INVALID_UINT),
      m_dataTransferId(LSF_INVALID_UINT) {
  }
  
  /// Standard c'tor.  Takes input values for all fields
  RunInfo( enums::Lsf::Platform p, enums::Lsf::DataOrigin d, UInt_t id, UInt_t startTime, UInt_t dataTransferId=LSF_INVALID_UINT )
    : m_platform(p), m_origin(d), m_id(id), m_startTime(startTime),
      m_dataTransferId(dataTransferId) {
  }
  
  /// Copy c'tor.  Nothing fancy, just copy all values.
  RunInfo( const RunInfo& other )
    : TObject(other), 
      m_platform(other.platform()), m_origin(other.dataOrigin()), 
      m_id(other.id()), m_startTime(other.startTime()),
      m_dataTransferId(other.dataTransferId()) {
  }
  
  /// D'tor.  Nothing special.
  virtual ~RunInfo() {
  }
  
  /// Assignement operator.  Nothing fancy, just copy all values.
  inline RunInfo& operator=( const RunInfo& other ) {
    initialize(other.platform(),other.dataOrigin(),
	       other.id(),other.startTime());
    setDataTransferId(other.dataTransferId());
    return *this;
  }
  
  /// The platform type this run was taken on.
  inline enums::Lsf::Platform platform() const { 
    return m_platform;
  }
  
  /// The type of data from this run (Orbit, MC or ground)
  inline enums::Lsf::DataOrigin dataOrigin() const {
    return m_origin;
  }
    
  /// The ground based ID of this run
  /// This is usually defined on the ground, but if the LAT DAQ reboots on-orbit, 
  /// The Ground ID can be reset
  inline UInt_t id() const {
    return m_id;
  }
  
  /// The start time of this run
  /// This is the number of seconds since GLAST epoch start when the run started
  /// FIXME (is this latched when the trigger was enabled?)
  inline UInt_t startTime() const {
    return m_startTime;
  }

  inline UInt_t dataTransferId() const {
     return m_dataTransferId;
  }
  
  /// set everything at once
  inline void initialize(enums::Lsf::Platform p, enums::Lsf::DataOrigin d, 
			 UInt_t id, UInt_t startTime,
                         UInt_t dataTransferId = LSF_INVALID_UINT ) {
    m_platform = p;
    m_origin = d;
    m_id = id;
    m_startTime = startTime;
    m_dataTransferId = dataTransferId;
  }
  
  // set the individual data members
  inline void setPlatform(enums::Lsf::Platform val) { m_platform = val; }
  inline void setDataOrigin(enums::Lsf::DataOrigin val) { m_origin = val; }
  inline void setId ( UInt_t val ) { m_id = val; }
  inline void setStartTime ( UInt_t val ) { m_startTime = val; }
  inline void setDataTransferId( UInt_t val ) { m_dataTransferId = val; }
  
  /// this is a helper function for formated printing
  void printString(std::ostream& s) const {
    static const char* Platforms[3] = {"Lat","Testbed","Host"};
    static const char* DataOrigins[3] = {"Orbit","MC","Ground"};
    if ( m_platform == enums::Lsf::NoPlatform ) {
      s << "None_";
    } else {
      s << Platforms[m_platform] << '_' ;
    }
    if ( m_origin == enums::Lsf::NoOrigin ) {
      s << "None_";
    } else {
      s << DataOrigins[m_origin] << '_';     
    }
    s << m_id << '_' << m_startTime;
    s << "dataTransferId: " << m_dataTransferId;
  }

  /// Reset function
  void Clear(Option_t* /* option="" */) {
    m_platform = enums::Lsf::NoPlatform;
    m_origin = enums::Lsf::NoOrigin;
    m_id = LSF_INVALID_UINT;
    m_startTime = LSF_INVALID_UINT;
    m_dataTransferId = LSF_INVALID_UINT;
  }

  /// ROOT print function
  void Print(Option_t* /* option="" */) const;

  /// Create a fake for tests
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ); 

  /// Compare to another in tests
  Bool_t CompareInRange( const RunInfo&, const std::string & name = "" ) const ; 

  /// For sorting, override TObject::Compare()
  virtual Int_t Compare(const TObject *obj) const; 
  
  /// Specify that this class is sortable
  virtual Bool_t IsSortable() const { return kTRUE; };


private:
  
  /// The platform type this run was taken on.
  enums::Lsf::Platform m_platform;
  
  /// The type of data from this run (Orbit, MC or ground)
  enums::Lsf::DataOrigin m_origin;
  
  /// The ground ID for the run
  UInt_t m_id;
  
  /// The time the run started (in seconds since the start of the GLAST epoch)
  UInt_t m_startTime;

  UInt_t m_dataTransferId;

  ClassDef(RunInfo,2) // Unique Run Identifier

};


#endif    // EVENT_RUNINFO_H

#ifndef LSF_RUNINFO_H
#define LSF_RUNINFO_H 1

#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "Event/TopLevel/Definitions.h"

#include "enums/Lsf.h"
/** @class RunInfo
* @brief encapsulate the run id parts of the event context
*
* 
*
* $Header$
*/

namespace LsfEvent {
    
  class RunInfo {
       
  public:

    /// Default c'tor.  Assigns sentinel values to all fields
    RunInfo( )
      : m_platform(enums::Lsf::NoPlatform), m_origin(enums::Lsf::NoOrigin), 
	m_id(LSF_INVALID_UINT), m_startTime(LSF_INVALID_UINT) {
    }
  
    /// Standard c'tor.  Takes input values for all fields
    RunInfo( enums::Lsf::Platform p, enums::Lsf::DataOrigin d, unsigned int id, unsigned int startTime )
      : m_platform(p), m_origin(d), m_id(id), m_startTime(startTime) {
    }
    
    /// Copy c'tor.  Nothing fancy, just copy all values.
    RunInfo( const RunInfo& other )
      : m_platform(other.platform()), m_origin(other.dataOrigin()), 
	m_id(other.id()), m_startTime(other.startTime()) {
    }
    
    /// D'tor.  Nothing special.
    virtual ~RunInfo() {
    }
    
    /// Assignement operator.  Nothing fancy, just copy all values.
    inline RunInfo& operator=( const RunInfo& other ) {
      set(other.platform(),other.dataOrigin(),
	  other.id(),other.startTime());
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
    inline unsigned int id() const {
      return m_id;
    }
    
    /// The start time of this run
    /// This is the number of seconds since GLAST epoch start when the run started
    /// FIXME (is this latched when the trigger was enabled?)
    inline unsigned int startTime() const {
      return m_startTime;
    }
    
    /// set everything at once
    inline void set(enums::Lsf::Platform p, enums::Lsf::DataOrigin d, 
		    unsigned int id, unsigned int startTime ) {
      m_platform = p;
      m_origin = d;
      m_id = id;
      m_startTime = startTime;
    }
    
    // set the individual data members
    inline void setPlatform(enums::Lsf::Platform val) { m_platform = val; }
    inline void setDataOrigin(enums::Lsf::DataOrigin val) { m_origin = val; }
    inline void setId ( unsigned int val ) { m_id = val; }
    inline void setStartTime ( unsigned int val ) { m_startTime = val; }
    
    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const RunInfo& obj )    {
      return s << obj.m_platform << ' ' << obj.m_origin << ' ' << obj.m_id << ' ' << obj.m_startTime;
    }
    /// Serialize the object for reading
    friend StreamBuffer& operator>> ( StreamBuffer& s, RunInfo& obj )          {
      return s << obj.m_platform << obj.m_origin << obj.m_id << obj.m_startTime;
    }
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const RunInfo& obj )    {
      return obj.fillStream(s);
    }

    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const                            {
      s << "class RunInfo : "
	<< EventField( EventFormat::field12 );
      printString(s);
      return s;
    }
    
  private:

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
    }
    
    /// The platform type this run was taken on.
    enums::Lsf::Platform m_platform;
    
    /// The type of data from this run (Orbit, MC or ground)
    enums::Lsf::DataOrigin m_origin;
    
    /// The ground ID for the run
    unsigned int m_id;
    
    /// The time the run started (in seconds since the start of the GLAST epoch)
    unsigned int m_startTime;
    
  };

}


#endif    // EVENT_RUNINFO_H

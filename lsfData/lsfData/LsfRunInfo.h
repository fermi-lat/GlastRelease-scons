#ifndef LSFDATA_RUNINFO_H
#define LSFDATA_RUNINFO_H 1

#include <iostream>
#include <iomanip>
#include <cstring>

#include "enums/Lsf.h"
/** @class RunInfo
* @brief encapsulate the run id parts of the event context
*
* 
*
* $Header$
*/

namespace lsfData {
    
  class RunInfo {
       
  public:
    
    RunInfo( )
      : m_platform(enums::Lsf::NoPlatform), m_origin(enums::Lsf::NoOrigin), m_id(0), m_startTime(0), m_dataTransferId(LSF_INVALID_UINT) {
    }
    
    RunInfo( enums::Lsf::Platform p, enums::Lsf::DataOrigin d, unsigned int id, unsigned int startTime, unsigned int dataTransferId = LSF_INVALID_UINT )
      : m_platform(p), m_origin(d), m_id(id), m_startTime(startTime),
        m_dataTransferId(dataTransferId) {
    }
    
    RunInfo( const RunInfo& other )
      : m_platform(other.platform()), m_origin(other.dataOrigin()), 
	m_id(other.id()), m_startTime(other.startTime()),
        m_dataTransferId(other.dataTransferId()) {
    }

    ~RunInfo() {
    }

    void clear() {
        m_platform = enums::Lsf::NoPlatform;
        m_origin = enums::Lsf::NoOrigin;
        m_id = 0;
        m_startTime = 0; 
        m_dataTransferId = LSF_INVALID_UINT;
    }

    void print() const {
      printf(" run:      groundid = 0x%08X, started = 0x%08X (%u)\n", id(), startTime(), startTime());
      printf(" run:      platform = (%d)\n",  platform());
      printf(" run:      origin = (%d)\n",  dataOrigin());
      printf(" dataTransferId: %u\n", dataTransferId());
    }

   /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const RunInfo& obj )    {
      return obj.fillStream(s);
    }

    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const {
      s << " run:      groundid = 0x" << std::hex << std::uppercase
        << std::setfill('0') << std::setw(8)
        << id() << std::dec << ", started = 0x" << std::hex 
        << startTime() << std::dec << " (" << startTime() << ")\n"
        << " run:      platform = (" << platform() << ")\n"
        << " run:      origin = (" << dataOrigin() << ")\n"
        << " dataTransferId: " << dataTransferId() << "\n";
      return s;
    }

    /// Assignement operator
    inline RunInfo& operator=( const RunInfo& other ) {
      set(other.platform(),other.dataOrigin(),
	  other.id(),other.startTime(),other.dataTransferId());
      return *this;
    }

    /// The platform this run was taken on
    inline enums::Lsf::Platform platform() const { 
      return m_platform;
    }
    
    /// The type of data from this run (Orbit, MC or ground) data
    inline enums::Lsf::DataOrigin dataOrigin() const {
      return m_origin;
    }
    
    /// The ground based ID of this run
    /// This is usually defined on the ground, but if the LAT DAQ reboots on-orbit, the 
    /// Ground ID can be reset
    inline unsigned int id() const {
      return m_id;
    }
    
    /// The start time of this run
    /// This is the number of seconds since GLAST epoch start when the run started
    /// FIXME (is this latched when the trigger was enabled?)
    inline unsigned int startTime() const {
      return m_startTime;
    }

    inline unsigned int dataTransferId() const {
        return m_dataTransferId; 
    }
    
    /// set everything at once
    inline void set(enums::Lsf::Platform p, enums::Lsf::DataOrigin d, 
		    unsigned int id, unsigned int startTime,
                    unsigned int dataTransferId=LSF_INVALID_UINT ) {
      m_platform = p;
      m_origin = d;
      m_id = id;
      m_startTime = startTime;
      m_dataTransferId = dataTransferId;
    }

    // set the individual data members
    inline void setPlatform(enums::Lsf::Platform val) { m_platform = val; }
    inline void setDataOrigin(enums::Lsf::DataOrigin val) { m_origin = val; }
    inline void setId ( unsigned int val ) { m_id = val; }
    inline void setStartTime ( unsigned int val ) { m_startTime = val; }
    inline void setDataTransferId( unsigned int val) { m_dataTransferId = val; }
    
    
  private:
    
    ///
    enums::Lsf::Platform m_platform;
    
    /// 
    enums::Lsf::DataOrigin m_origin;
    
    ///
    unsigned int m_id;
    
    /// storage for ctx.run.startedAt
    /// flight run id
    unsigned int m_startTime;

    /// created to store the LSFReader::runId obtained from evt.  This is
    /// set at start-of-decode time from a value provided in the RetDef
    /// file for each decoding segment, so it varies depending on from which
    /// downlink the event was decoded.
    unsigned int m_dataTransferId;
  };

}


#endif    // EVENT_RUNINFO_H

#ifndef LSF_DATAGRAMINFO_H
#define LSF_DATAGRAMINFO_H 1

#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "Event/TopLevel/Definitions.h"

#include "enums/Lsf.h"

/** @class DatagramInfo
* @brief encapsulate the Datagram id parts of the event context
*
* 
*
* $Header$
*/

namespace LsfEvent {
    
  class DatagramInfo {
    
  public:

    DatagramInfo( )
      : m_openReason(enums::Lsf::NoOpenReason), m_openRequester(enums::Lsf::NoOpenRequest), 
	m_crate(enums::Lsf::NoCrate), m_mode(enums::Lsf::NoMode), 
	m_closeReason(enums::Lsf::NoCloseReason), m_closeRequester(enums::Lsf::NoCloseRequest),
	m_datagrams(0),m_modeChanges(0){
    };
        
    DatagramInfo( enums::Lsf::OpenReason oReason, enums::Lsf::OpenRequester oReq, 
		  enums::Lsf::Crate c, enums::Lsf::Mode m, 
		  enums::Lsf::CloseReason cReason, enums::Lsf::CloseRequester cReq,
		  unsigned int datagrams, unsigned int modeChanges)
      : m_openReason(oReason), m_openRequester(oReq), 
	m_crate(c), m_mode(m), 
	m_closeReason(cReason), m_closeRequester(cReq),
	m_datagrams(datagrams),m_modeChanges(modeChanges){
    }
    
    ~DatagramInfo() {
    }

    DatagramInfo( const DatagramInfo& other ) 
      : m_openReason(other.openReason()), m_openRequester(other.openRequester()),
	m_crate(other.crate()), m_mode(other.mode()), 
	m_closeReason(other.closeReason()), m_closeRequester(other.closeRequester()),
	m_datagrams(other.datagrams()),m_modeChanges(other.modeChanges()){
    }
    
    inline DatagramInfo& operator=(const DatagramInfo& other) {
      set(other.openReason(),other.openRequester(),
	  other.crate(),other.mode(),
	  other.closeReason(),other.closeRequester(),
	  other.datagrams(),other.modeChanges());
      return *this;
    }

    /// number of mode changes since the start of the run
    inline unsigned int modeChanges() const { return m_modeChanges; };

    /// number of datagrams sent since the start of the run
    /// this is identical to the datagram sequence number
    inline unsigned int datagrams() const { return m_datagrams; };

    /// Reason that this datagram was opened
    /// FIXME: 
    inline enums::Lsf::OpenReason openReason() const { return m_openReason; }

    /// Method used to open this datagram
    inline enums::Lsf::OpenRequester openRequester() const { return m_openRequester; }

    /// Source that this datagram came from
    inline enums::Lsf::Crate crate() const { return m_crate; }

    /// Operating mode the LAT was in when the data for this data were acquired
    inline enums::Lsf::Mode mode() const { return m_mode; } 

    /// Reason that this datagram was opened
    /// FIXME: should this be LsfEvent::DatagramInfo::Stop for the last datagram in a run?
    inline enums::Lsf::CloseReason closeReason() const { return m_closeReason; }

    /// Method used to close this datagram
    inline enums::Lsf::CloseRequester closeRequester() const { return m_closeRequester; }

    /// set everything at once
    inline void set(enums::Lsf::OpenReason oReason, enums::Lsf::OpenRequester oReq, 
		    enums::Lsf::Crate c, enums::Lsf::Mode m, 
		    enums::Lsf::CloseReason cReason, enums::Lsf::CloseRequester cReq,
		    unsigned int datagrams, unsigned int modeChanges) {
      m_openReason = oReason;
      m_openRequester = oReq;
      m_crate = c;
      m_mode = m;
      m_closeReason = cReason;
      m_closeRequester = cReq;
      m_datagrams = datagrams;
      m_modeChanges = modeChanges;
    }
    
    // set the individual data members
    inline void setModeChanges( unsigned int val ) { m_modeChanges = val; };
    inline void setDatagrams( unsigned int val ) { m_datagrams = val; }; 
    inline void setOpenReason( enums::Lsf::OpenReason val ) { m_openReason = val; };
    inline void setOpenRequester( enums::Lsf::OpenRequester val ) { m_openRequester = val; };
    inline void setCrate( enums::Lsf::Crate val ) { m_crate = val; };
    inline void setMode( enums::Lsf::Mode val ) { m_mode = val; }; 
    inline void setCloseReason( enums::Lsf::CloseReason val ) { m_closeReason = val; };
    inline void setCloseRequester( enums::Lsf::CloseRequester val ) { m_closeRequester = val; };

    
    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const DatagramInfo& obj )    {
      return s << obj.m_openReason << ' ' << obj.m_openRequester << ' ' 
	       << obj.m_crate << ' ' << obj.m_mode << ' '
	       << obj.m_closeReason << ' ' << obj.m_closeRequester << ' ' 
	       << obj.m_modeChanges << ' ' << obj.m_datagrams;
	}
    /// Serialize the object for reading
    friend StreamBuffer& operator>> ( StreamBuffer& s, DatagramInfo& obj )          {
      // FIXME StreamBuffer can't read enums
      //return s >> obj.m_openReason  >> obj.m_openRequester 
      //       >> obj.m_crate >> obj.m_mode 
      //       >> obj.m_closeReason >> obj.m_closeRequester
      //       >> obj.m_modeChanges>> obj.m_datagrams;
      return s;
    }
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const DatagramInfo& obj )    {
      return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const                            {
      s << "class DatagramInfo : "
	<< EventField( EventFormat::field12 )
	<< m_openReason << ' ' << m_openRequester << ' ' 
	<< m_crate << ' ' << m_mode << ' '
	<< m_closeReason << ' ' << m_closeRequester << ' ' 
	<< m_modeChanges << ' ' << m_datagrams;      
      return s;
    }
    
  private:

    enums::Lsf::OpenReason m_openReason;
    enums::Lsf::OpenRequester m_openRequester;
    enums::Lsf::Crate m_crate;
    enums::Lsf::Mode m_mode;    
    enums::Lsf::CloseReason m_closeReason;
    enums::Lsf::CloseRequester m_closeRequester;

    unsigned int m_datagrams;    
    unsigned int m_modeChanges; 

  };

}


#endif    // EVENT_DATAGRAMINFO_H

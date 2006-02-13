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
    /// Default c'tor.  Assigns sentinel values to all fields
    DatagramInfo( )
      : m_openAction(enums::Lsf::Open::Unspecified), m_openReason(enums::Lsf::Open::Unknown), 
	m_crate(enums::Lsf::NoCrate), m_mode(enums::Lsf::NoMode), 
	m_closeAction(enums::Lsf::Close::Unspecified), m_closeReason(enums::Lsf::Close::Unknown), 
	m_datagrams(LSF_INVALID_UINT),m_modeChanges(LSF_INVALID_UINT){
    };
    
    /// Standard c'tor.  Takes input values for all fields
    DatagramInfo( enums::Lsf::Open::Action oAction, enums::Lsf::Open::Reason oReason, 
		  enums::Lsf::Crate c, enums::Lsf::Mode m, 
		  enums::Lsf::Close::Action cAction, enums::Lsf::Close::Reason cReason, 
		  unsigned int datagrams, unsigned int modeChanges)
      : m_openAction(oAction), m_openReason(oReason), 
	m_crate(c), m_mode(m), 
	m_closeAction(cAction), m_closeReason(cReason),
	m_datagrams(datagrams),m_modeChanges(modeChanges){
    }

    /// Copy c'tor.  Nothing fancy, just copy all values.
    DatagramInfo( const DatagramInfo& other ) 
      : m_openAction(other.openAction()), m_openReason(other.openReason()),
	m_crate(other.crate()), m_mode(other.mode()), 
	m_closeAction(other.closeAction()), m_closeReason(other.closeReason()),
	m_datagrams(other.datagrams()),m_modeChanges(other.modeChanges()){
    }
    
    /// D'tor.  Nothing special.
    virtual ~DatagramInfo() {
    }
    
    /// Assignement operator.  Nothing fancy, just copy all values.
    inline DatagramInfo& operator=(const DatagramInfo& other) {
      set(other.openAction(),other.openReason(),
	  other.crate(),other.mode(),
	  other.closeAction(),other.closeReason(),
	  other.datagrams(),other.modeChanges());
      return *this;
    }
    
    /// number of mode changes since the start of the run
    inline unsigned int modeChanges() const { return m_modeChanges; };
    
    /// number of datagrams sent since the start of the run (aka the datagram sequence number)
    inline unsigned int datagrams() const { return m_datagrams; };
    
    /// The action that caused the datagram to be opened
    inline enums::Lsf::Open::Action openAction() const { return m_openAction; }
    
    /// The reason this datagram was opened
    inline enums::Lsf::Open::Reason openReason() const { return m_openReason; }
    
    /// Source that this datagram came from
    inline enums::Lsf::Crate crate() const { return m_crate; }
    
    /// Operating mode the LAT was in when the data for this datagram were acquired
    inline enums::Lsf::Mode mode() const { return m_mode; } 
    
    /// The action that caused the datagram to be closed
    inline enums::Lsf::Close::Action closeAction() const { return m_closeAction; }
    
    /// The reason this datagram was closed
    inline enums::Lsf::Close::Reason closeReason() const { return m_closeReason; }
    
    /// set everything at once
    inline void set(enums::Lsf::Open::Action oAction, enums::Lsf::Open::Reason oReason,
		    enums::Lsf::Crate c, enums::Lsf::Mode m, 
		    enums::Lsf::Close::Action cAction, enums::Lsf::Close::Reason cReason,
		    unsigned int datagrams, unsigned int modeChanges) {
      m_openAction = oAction;
      m_openReason = oReason;
      m_crate = c;
      m_mode = m;
      m_closeAction = cAction;
      m_closeReason = cReason;
      m_datagrams = datagrams;
      m_modeChanges = modeChanges;
    }
    
    // set the individual data members
    inline void setModeChanges( unsigned int val ) { m_modeChanges = val; };
    inline void setDatagrams( unsigned int val ) { m_datagrams = val; }; 
    inline void setOpenAction( enums::Lsf::Open::Action val ) { m_openAction = val; };
    inline void setOpenReason( enums::Lsf::Open::Reason val ) { m_openReason = val; };
    inline void setCrate( enums::Lsf::Crate val ) { m_crate = val; };
    inline void setMode( enums::Lsf::Mode val ) { m_mode = val; }; 
    inline void setCloseAction( enums::Lsf::Close::Action val ) { m_closeAction = val; };
    inline void setCloseReason( enums::Lsf::Close::Reason val ) { m_closeReason = val; };
    
    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const DatagramInfo& obj )    {
      return s << obj.m_openAction << ' ' << obj.m_openReason << ' ' 
	       << obj.m_crate << ' ' << obj.m_mode << ' '
	       << obj.m_closeAction << ' ' << obj.m_closeReason << ' '
	       << obj.m_modeChanges << ' ' << obj.m_datagrams;
	}
    /// Serialize the object for reading
    friend StreamBuffer& operator>> ( StreamBuffer& s, DatagramInfo& obj )          {
      // FIXME StreamBuffer can't read enums
      //return s >>  obj.m_openAction >> obj.m_openReason 
      //       >> obj.m_crate >> obj.m_mode 
      //       >> obj.m_closeAction >> obj.m_closeReason
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
	<< m_openAction << ' ' << m_openReason << ' ' 
	<< m_crate << ' ' << m_mode << ' '
	<< m_closeAction << ' ' << m_closeReason << ' ' 
	<< m_modeChanges << ' ' << m_datagrams;      
      return s;
    }
    
  private:

    /// The action that caused the datagram to be opened
    enums::Lsf::Open::Action m_openAction;
    /// The reason this datagram was opened
    enums::Lsf::Open::Reason m_openReason;
    
    /// Source that this datagram came from
    enums::Lsf::Crate m_crate;
    /// Operating mode the LAT was in when the data for this datagram were acquired
    enums::Lsf::Mode m_mode;    
    
    /// The action that caused the datagram to be closed
    enums::Lsf::Close::Action m_closeAction;
    /// The reason this datagram was closed
    enums::Lsf::Close::Reason m_closeReason;
    
    /// number of datagrams sent since the start of the run (aka the datagram sequence number)
    unsigned int m_datagrams;    
    /// number of mode changes since the start of the run
    unsigned int m_modeChanges; 
    
  };

}


#endif    // EVENT_DATAGRAMINFO_H

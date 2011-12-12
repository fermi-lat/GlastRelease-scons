#ifndef LHCBEVENT_MCEVENT_H
#define LHCBEVENT_MCEVENT_H 1

#include <iostream>
#include <vector>
#include <algorithm>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/StreamBuffer.h"
//#include "GlastEvent/TopLevel/SubMCEvent.h"
#include "Event/TopLevel/Definitions.h"
#include "Event/Utilities/TimeStamp.h"
#include "GaudiKernel/IInterface.h"

static const CLID& CLID_McEvent = InterfaceID("McEvent", 1, 1);

namespace Event {

/** @class MCEvent
* @brief Top level Monte Carlo data object
* It can be identified by "/Event/MC" on the TDS
* 
* It contains:
* - run number
* - sequence number for keying random number generator
* - source ID 
* 
* $Header$
*/

class MCEvent : public DataObject                                              {
    
public:
    MCEvent( ):  m_sourceId(-1), m_run(0), m_sequence(0), m_sourceName("_") {}
    
    virtual ~MCEvent()  { }
    
    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const { return MCEvent::classID(); }
    static const CLID& classID() { return CLID_McEvent; }
    
    /// Clone operator
    MCEvent& operator=(const MCEvent& copy)                                      {
        return *this;
    }
    
    /// Retrieve 
    int getSourceId() const { return m_sourceId; }
    int getRunNumber() const      { return m_run; }
    int getSequence() const { return m_sequence; }
    const std::string& getSourceName()const{return m_sourceName;}

    /// initialize
    void initialize(int run, int source, long int seq, TimeStamp time, std::string name="_") {
        m_run = run; m_sourceId = source; m_sequence = seq; m_time=time; m_sourceName = name;}

    void setSourceName(const std::string& name){m_sourceName=name;}

    /// Retrieve reference to event time stamp
    const TimeStamp& time () const                              { return m_time; }
    /// Update reference to event time stamp
    void setTime (const TimeStamp& value)                      { m_time = value; }
    
    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const MCEvent& obj )      {
        return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    virtual std::ostream& fillStream( std::ostream& s ) const;
    
private:
    // identifier of the source
    int m_sourceId;
    /// run number
    unsigned int m_run;
    /// sequence number
    unsigned int m_sequence;
    
    /// Time stamp: use special class to encapsulate type
    TimeStamp           m_time;

    std::string m_sourceName; ///< "name for the source, hopefully unique"

};

//
// Inline code must be outside the class definition


/// Serialize the object for writing
inline StreamBuffer& MCEvent::serialize( StreamBuffer& s ) const               {
    DataObject::serialize(s);
    s << m_sourceId << m_sequence << m_run << m_time;
    return s;
}


/// Serialize the object for reading
inline StreamBuffer& MCEvent::serialize( StreamBuffer& s )                     {
    DataObject::serialize(s);
    s >> m_time >> m_run >> m_sequence >> m_sourceId;
    
    return s;
}


/// Fill the ASCII output stream
inline std::ostream& MCEvent::fillStream( std::ostream& s ) const              {
    s << "class MCEvent :\n"
        << "    Source Id = "
        << EventField( EventFormat::field12 ) << m_sourceId
        << "    Run number = "
        << EventField( EventFormat::field12 )  << m_run
        << "    Sequence = "
        << EventField( EventFormat::field12 )  << m_sequence
        << "    TimeStamp = "
        << EventField( EventFormat::field12 )  << m_time
        ;
    return s;
}

} // namespace Event
#endif    // GLASTEVENT_MCEVENT_H

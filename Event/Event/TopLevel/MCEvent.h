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

extern const CLID& CLID_McEvent;
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
    MCEvent( ) {}
    
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

    /// initialize
    void initialize(int run, int source, long int seq) {
        m_run = run; m_sourceId = source; m_sequence = seq;}
    
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
    
};

//
// Inline code must be outside the class definition


/// Serialize the object for writing
inline StreamBuffer& MCEvent::serialize( StreamBuffer& s ) const               {
    DataObject::serialize(s);
    s << m_sourceId << m_sequence << m_run;
    return s;
}


/// Serialize the object for reading
inline StreamBuffer& MCEvent::serialize( StreamBuffer& s )                     {
    DataObject::serialize(s);
    s >> m_run >> m_sequence >> m_sourceId;
    
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
        ;
    return s;
}

} // namespace Event
#endif    // GLASTEVENT_MCEVENT_H

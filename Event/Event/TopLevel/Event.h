#ifndef GLASTEVENT_EVENT_H
#define GLASTEVENT_EVENT_H 1

#include <iostream>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/DataObject.h"
#include "GlastEvent/Utilities/TimeStamp.h"
#include "GlastEvent/TopLevel/Definitions.h"

extern const CLID& CLID_Event;

/** @class Event
* @brief Essential header information of the event.
* It can be identified by "/Event" on the TDS.
*
* It contains:
* - run number
* - event number
* - time stamp
*
* $Header$
*/

class Event : public DataObject                                                {
    
public:
    
    Event()
        : DataObject() { }
    
    virtual ~Event() { }
    
    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const { return Event::classID(); }
    static const CLID& classID() { return CLID_Event; }
    
    /// Retrieve event number
    long event () const                                        { return m_event; }
    /// Update event number
    void setEvent (long value)                                { m_event = value; }
    
    /// Retrieve run number
    long run () const                                            { return m_run; }
    /// Update run number
    void setRun (long value)                                    { m_run = value; }
    
    /// Retrieve reference to event time stamp
    const TimeStamp& time () const                              { return m_time; }
    /// Update reference to event time stamp
    void setTime (const TimeStamp& value)                      { m_time = value; }
    
    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const Event& obj )        {
        return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    virtual std::ostream& fillStream( std::ostream& s ) const;
    
private:
    /// Event number
    long                m_event;
    /// Run number
    long                m_run;
    /// Time stamp
    TimeStamp           m_time;
};


//
// Inline code must be outside the class definition
//


/// Serialize the object for writing
inline StreamBuffer& Event::serialize( StreamBuffer& s ) const                 {
    DataObject::serialize(s);
    // HMK The member variables are not filled yet
    return s;
    //<< m_event
    //<< m_run
    //<< m_time;
}


/// Serialize the object for reading
inline StreamBuffer& Event::serialize( StreamBuffer& s )                       {
    DataObject::serialize(s);
    // HMK The member variables are not filled yet
    return s;
    //>> m_event
    //>> m_run
    //>> m_time;
}


/// Fill the output stream (ASCII)
inline std::ostream& Event::fillStream( std::ostream& s ) const                {
    return s
        << "class Event :"
        << "\n    Event number = "
        << GlastEventField( GlastEvent::field12 )
        << m_event
        << "\n    Run number   = "
        << GlastEventField( GlastEvent::field12 )
        << m_run
        << "\n    Time         = " << m_time;
}


#endif    // GLASTEVENT_EVENT_H

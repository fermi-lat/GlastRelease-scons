#ifndef Event_EVENTOVERLAY_H
#define Event_EVENTOVERLAY_H 1

#include <iostream>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"
#include "Event/Utilities/TimeStamp.h"
#include "Event/TopLevel/Definitions.h"
#include "enums/EventFlags.h"

static const CLID& CLID_EventOverlay = InterfaceID("EventOverlay", 1, 0);

/** @class EventOverlay
* @brief Essential header information of the event.
* It can be identified by "/Event" on the TDS.
*
* It contains:
* - run number
* - event number
* - time stamp
* - trigger word
* - trigger engines
* - trigger prescales
* - flag if the prescale expired
*
* $Header$
*/
namespace Event{

class EventOverlay : public DataObject                                                {
    
public:
    /// default constructor: flag all fields invalic
    EventOverlay()
        : DataObject(), m_event(0), m_run(-1), m_time(-1),  m_trigger((unsigned int)-1),
        m_triggerWordTwo((unsigned int)-1), m_gemPrescale((unsigned int)-1), 
        m_gltPrescale((unsigned int)-1), m_prescaleExpired((unsigned int)-1),
        m_gleamEventFlags(0) { }
    
    virtual ~EventOverlay() { }
    
    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const { return EventOverlay::classID(); }
    static const CLID& classID() { return CLID_EventOverlay; }
    
    /// Retrieve event number
    unsigned long long event () const        { return m_event; }
    /// Update event number
    void setEvent (unsigned long long value) { m_event = value; }
    
    /// Retrieve run number
    long run () const           { return m_run; }
    /// Update run number
    void setRun (long value)    { m_run = value; }
    
    /// Retrieve reference to event time stamp
    const TimeStamp& time () const                { return m_time; }
    /// Update reference to event time stamp
    void setTime (const TimeStamp& value)         { m_time = value; }

    /// Retrieve trigger word
    unsigned int trigger()const         { return m_trigger;}
    /// update trigger word
    void setTrigger(unsigned int value) {m_trigger = value;}

    /// Retrieve second trigger word
    unsigned int triggerWordTwo() const             {return m_triggerWordTwo; }
    /// update second trigger word
    void setTriggerWordTwo(unsigned int value)      { m_triggerWordTwo = value; }
    
    /// Retrieve the GEM prescale
    unsigned int gemPrescale() const        {return m_gemPrescale; }
    /// update the GEM prescale
    void setGemPrescale(unsigned int value) { m_gemPrescale = value; }
 
    /// Retrieve the GLT prescale
    unsigned int gltPrescale() const        {return m_gltPrescale; }
    /// update second trigger word
    void setGltPrescale(unsigned int value) { m_gltPrescale = value; }
 
    /// Retrieve prescale expired flag
    unsigned int prescaleExpired() const        {return m_prescaleExpired; }
    /// update second trigger word
    void setPrescaleExpired(unsigned int value) { m_prescaleExpired = value; }
    
    /// Retrieve live time
    double livetime()const         { return m_livetime;}
    /// update live time
    void setLivetime(double value) {m_livetime = value;}


    unsigned int gleamEventFlags() const    { return m_gleamEventFlags; }
    void setGleamEventFlags(unsigned int f) { m_gleamEventFlags = f; }
    void setTkrReconError()                 { m_gleamEventFlags |= enums::Gleam::TKRRECON; }
    void setCalReconError()                 { m_gleamEventFlags |= enums::Gleam::CALRECON; }
    void setAcdReconError()                 { m_gleamEventFlags |= enums::Gleam::ACDRECON; }
    void setAnalysisNtupleError()           { m_gleamEventFlags |= enums::Gleam::ANALYSISNTUPLE; }
    bool goodEvent()                        { return (m_gleamEventFlags == 0); }
    bool badEvent()                         { return (m_gleamEventFlags != 0); }
    bool badTkrRecon()                      { return (m_gleamEventFlags & enums::Gleam::TKRRECON)       != 0; }
    bool badAcdRecon()                      { return (m_gleamEventFlags & enums::Gleam::ACDRECON)       != 0; }
    bool badCalRecon()                      { return (m_gleamEventFlags & enums::Gleam::CALRECON)       != 0; }
    bool badAnalysisNtuple()                { return (m_gleamEventFlags & enums::Gleam::ANALYSISNTUPLE) != 0; }
 
    
    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const EventOverlay& obj )        {
        return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    virtual std::ostream& fillStream( std::ostream& s ) const;
    
private:
    /// Event number
    unsigned long long  m_event;
    /// Run number
    long                m_run;
    /// Time stamp: use special class to encapsulate type
    TimeStamp           m_time;
    /// trigger word: note that a class is available to use for this (THB: do we need this?)
    unsigned int        m_trigger;
    /// store additional bits, now that we're saving the trigger engines for both GEM and GLT
    unsigned int        m_triggerWordTwo;

    /// live time
    double              m_livetime;
    /// GEM trigger prescale
    unsigned int        m_gemPrescale;
    /// GLT trigger prescale
    unsigned int        m_gltPrescale;
    /// flag if prescale counter expired on this event
    unsigned int        m_prescaleExpired;

    /// error bits, where m_gleamEventFlags == 0 is a GOOD event
    unsigned int        m_gleamEventFlags;
};


//
// Inline code must be outside the class definition
//


/// Serialize the object for writing
inline StreamBuffer& EventOverlay::serialize( StreamBuffer& s ) const                 {
    DataObject::serialize(s);
    // HMK The member variables are not filled yet
    return s;
    //<< m_event
    //<< m_run
    //<< m_time;
}


/// Serialize the object for reading
inline StreamBuffer& EventOverlay::serialize( StreamBuffer& s )                       {
    DataObject::serialize(s);
    // HMK The member variables are not filled yet
    return s;
    //>> m_event
    //>> m_run
    //>> m_time;
}


/// Fill the output stream (ASCII)
inline std::ostream& EventOverlay::fillStream( std::ostream& s ) const                {
    return s
        << "class Event :"
        << "\n    Event number = "
        << EventField( EventFormat::field12 )
        << m_event
        << "\n    Run number   = "
        << EventField( EventFormat::field12 )
        << m_run
        << "\n    Time         = " << m_time
        << std::setbase(16) << "\n    Trigger      = " << m_trigger
        << "\n    TriggerWordTwo      = " << m_triggerWordTwo
        << std::setbase(10) << "\n    GEM prescale = " << m_gemPrescale
        << "\n    GLT prescale = " << m_gltPrescale
        << "\n    Prescale expired    = " << m_prescaleExpired;
}
} // namespace Event
#endif    // Event_EVENT_H

#ifndef LHCBEVENT_MCEVENT_H
#define LHCBEVENT_MCEVENT_H 1

#include <iostream>
#include <vector>
#include <algorithm>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/SubMCEvent.h"
#include "GlastEvent/TopLevel/Definitions.h"

extern const CLID& CLID_McEvent;

/** @class MCEvent
* @brief Top level Monte Carlo data object
* It can be identified by "/Event/MC" on the TDS
* 
* It contains:
* - source ID
* 
* $Header$
*/

class MCEvent : public DataObject                                              {
    
public:
    MCEvent( int sourceId=0)
        : DataObject(), m_sourceId(sourceId) {}
    
    virtual ~MCEvent()  { }
    
    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const { return MCEvent::classID(); }
    static const CLID& classID() { return CLID_McEvent; }
    
    /// Clone operator
    MCEvent& operator=(const MCEvent& copy)                                      {
        return *this;
    }
    
    /// Retrieve 
    int sourceId () const                                                         {
        return m_sourceId;
    }
    /// set
    void setSourceId(int s);
    
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
    
};

inline void MCEvent::setSourceId(int s)
{
    m_sourceId = s;
}
//
// Inline code must be outside the class definition
//
//#include "GlastEvent/MonteCarlo/MCVertex.h"  HMA commented this out until we have MCVertex.h


/// Serialize the object for writing
inline StreamBuffer& MCEvent::serialize( StreamBuffer& s ) const               {
    DataObject::serialize(s);
    s << m_sourceId;
    return s;
}


/// Serialize the object for reading
inline StreamBuffer& MCEvent::serialize( StreamBuffer& s )                     {
    DataObject::serialize(s);
    s >> m_sourceId;
    return s;
}


/// Fill the ASCII output stream
inline std::ostream& MCEvent::fillStream( std::ostream& s ) const              {
    s << "class MCEvent :\n"
        << "    Source Id = "
        << GlastEventField( GlastEvent::field12 )
        << m_sourceId;
    return s;
}


#endif    // GLASTEVENT_MCEVENT_H

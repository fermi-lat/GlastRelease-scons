#ifndef GLASTEVENT_DIGIEVENT_H
#define GLASTEVENT_DIGIEVENT_H 1

#include <iostream>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/DataObject.h"
#include "GlastEvent/Utilities/TriggerPattern.h"
#include "GlastEvent/TopLevel/Definitions.h"

extern const CLID& CLID_DigiEvent;

/** @class DigiEvent
* @brief Defines the top level object for digitization data.
* It can be identified by "/Event/Digi" on the TDS
* 
* It contains:
* - flag, if comming from Monte Carlo
* 
*/

class DigiEvent : public DataObject                                             {
    
public:
    
    DigiEvent()
        : DataObject(), m_fromMC(false) { }
    
    virtual ~DigiEvent() { }
    
    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const  { return DigiEvent::classID(); }
    static const CLID& classID() { return CLID_DigiEvent; }
    
    ///  Retrieve flag of origin
    bool fromMC () const                                                         {
        return m_fromMC;
    }
    ///  Update flag of origine
    void setFromMC (bool value)                                                  {
        m_fromMC = value;
    }
    
    /// Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    /// Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const DigiEvent& obj )     {
        return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    virtual std::ostream& fillStream( std::ostream& s ) const;
    
private: 
    /// Flag of origin
    bool m_fromMC;
};


//
// Inline code must be outside the class definition
//


/// Serialize the object for writing
inline StreamBuffer& DigiEvent::serialize( StreamBuffer& s ) const              {
    DataObject::serialize(s);
    unsigned char u = (m_fromMC) ? 1 : 0;
    return s << u;
}


/// Serialize the object for reading
inline StreamBuffer& DigiEvent::serialize( StreamBuffer& s )                    {
    DataObject::serialize(s);
    unsigned char u;
    s >> u;
    m_fromMC = (u) ? true : false;
    return s;
}


/// Fill the output stream (ASCII)
inline std::ostream& DigiEvent::fillStream( std::ostream& s ) const             {
    s << "class DigiEvent :"
        << "\n    Flag of origin    = ";
    if( m_fromMC ) {
        s << " true";
    }
    else {
        s << "false";
    }
    return s;
}


#endif  // GLASTEVENT_DIGIEVENT_H

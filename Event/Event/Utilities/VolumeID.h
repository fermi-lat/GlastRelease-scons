// $Header$
#ifndef GlastEvent_VolumeID_H
#define GlastEvent_VolumeID_H 1

// If you wish to introduce the namespace `GlastEvent', uncomment
// the lines commented as `NameSpace'.


// Include files
#include <iostream>
#include "Gaudi/Kernel/StreamBuffer.h"
#include "GlastEvent/TopLevel/Definitions.h"

/*!
  Class: VolumeID
  Description: The base class for GLAST media volume identification
  Author: OZAKI Masanobu
  History:
    2000-12-07 M.Ozaki : created
 */

//namespace GlastEvent { // NameSpace

class VolumeID {
public:
    /// Constructor
    VolumeID(){};
    /// Destructor
    ~VolumeID(){};

    /// Retrieve volume ID
    long volumeID() const;
    /// Update volume ID
    void setVolumeID( long value );

    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const VolumeID& obj ){
        return s << obj.m_volumeID;
    }
    /// Serialize the object for reading
    friend StreamBuffer& operator>> ( StreamBuffer& s, VolumeID& obj ){
        return s >> obj.m_volumeID;
    }

    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const VolumeID& obj ){
        return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const{
        return s << "class VolumeID : "
          << GlastEventField( GlastEvent::field4 )
          << m_volumeID;
    }


private:
    long  m_volumeID;
};

//} // NameSpace GlastEvent


// Inline codes

//namespace GlastEvent { // NameSpace

/// Retrieve volume ID
inline long VolumeID::volumeID() const
{
    return m_volumeID;
}


/// Update volume ID
inline void VolumeID::setVolumeID( long value )
{
    m_volumeID = value;
}

//} // NameSpace GlastEvent

#endif // GlastEvent_VolumeID_H

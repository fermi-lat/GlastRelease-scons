// $Header$
#ifndef GlastEvent_VolumeID_H
#define GlastEvent_VolumeID_H 1


/*!
  Class: VolumeID
  Description: The base class for GLAST media volume identification
  Author: OZAKI Masanobu
  History:
    2000-12-07 M.Ozaki : created
 */


class VolumeID {
public:
    /// Constructor
    VolumeID();
    /// Destructor
    ~VolumeID();

    /// Retrieve volume ID
    long volumeID() const;
    /// Update volume ID
    long setVolumeID( long value );

private:
    long  m_volumeID;
};


#include "GlastEvent/Utilities/VolumeID.cpp"


#endif // GlastEvent_VolumeID_H

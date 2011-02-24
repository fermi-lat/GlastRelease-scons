/** @file ITkrTrackVecTool.h
@brief Abstract interface to return a vector of pointers to *all* the tracks
@author Leon Rochester
$Header$
*/

#ifndef _H_ITkrTrackVecTool
#define _H_ITkrTrackVecTool

#include "GaudiKernel/IAlgTool.h"

#include "Event/Recon/TkrRecon/TkrTrack.h"

#include <vector>

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_ITkrTrackVecTool("ITkrTrackVecTool", 0 , 0); 

/** @class ITkrTrackVecTool
* @brief Abstract interface: returns a vector of pointers to *all* the tracks
*/

class   ITkrTrackVecTool : virtual public IAlgTool {
public:

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ITkrTrackVecTool; }

    virtual std::vector<Event::TkrTrack*> getTrackVec() = 0;
};

#endif  // _H_ITkrTrackVecTool

#ifndef IAcdFailureModeSvc_H
#define IAcdFailureModeSvc_H

// $Header$

/** @file
@author H. Kelly 
*/


#include "GaudiKernel/IInterface.h"
#include "idents/AcdId.h"

// Declaration of the interface ID ( interface id, major version,

// minor version)

static const InterfaceID IID_IAcdFailureModeSvc("IAcdFailureModeSvc", 1 , 0);



/** @class IAcdFailureModeSvc
* @brief Interface class for AcdFailureModeSvc
*
* Author:  H. Kelly
*
*/

class IAcdFailureModeSvc : virtual public IInterface {

public:

    static const InterfaceID& interfaceID() { 
        return IID_IAcdFailureModeSvc; 
    }

    /// look for AcdId in list of dead tiles
    virtual bool matchAcdId(idents::AcdId id)=0;

protected:

    /// process the input list of Tile Ids
    virtual void processDetectorList()=0;

};

#endif // IAcdFailureModeSvc_H


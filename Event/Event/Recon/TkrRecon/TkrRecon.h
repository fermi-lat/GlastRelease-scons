#ifndef Event_TKRRECON_H
#define Event_TKRRECON_H 1

#include <iostream>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/DataObject.h"
#include "Event/TopLevel/Definitions.h"

#include "GaudiKernel/IInterface.h"

static const CLID& CLID_TkrRecon = InterfaceID("TkrRecon", 1, 0);

/** @class TkrRecon
* @brief Defines the top level object for digitization data.
* It can be identified by "/Event/TkrRecon" on the TDS
* 
* 
* $Header$
*/
namespace Event {  // NameSpace

class TkrRecon : public DataObject {
    
public:
    
    TkrRecon()
        : DataObject() { }
    
    virtual ~TkrRecon() { }
    

    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const  { return TkrRecon::classID(); }
    static const CLID& classID() { return CLID_TkrRecon; }
    
};

}


#endif  // GLASTEVENT_TKRRECON_H

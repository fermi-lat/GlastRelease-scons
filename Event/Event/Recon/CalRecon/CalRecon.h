#ifndef Event_CALRECON_H
#define Event_CALRECON_H 1

#include <iostream>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/DataObject.h"
#include "Event/TopLevel/Definitions.h"

extern const CLID& CLID_CalRecon;

/** @class CalRecon
* @brief Defines the top level object for digitization data.
* It can be identified by "/Event/CalRecon" on the TDS
* 
* 
* $Header$
*/
namespace Event {  // NameSpace

class CalRecon : public DataObject {
    
public:
    
    CalRecon()
        : DataObject() { }
    
    virtual ~CalRecon() { }
    

    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const  { return CalRecon::classID(); }
    static const CLID& classID() { return CLID_CalRecon; }
    
};

}


#endif  // GLASTEVENT_CALRECON_H

#ifndef ICalFailureModeSvc_H
#define ICalFailureModeSvc_H 1

// Include files
#include "GaudiKernel/IInterface.h"
#include "idents/CalXtalId.h"
#include <vector>
#include <map>

// Declaration of the interface ID ( interface id, major version,
// minor version)
static const InterfaceID IID_ICalFailureModeSvc("ICalFailureModeSvc", 1 , 0);

/** @class ICalFailureModeSvc
* @brief Interface class for CalFailureModeSvc
*
* Author:  R.Dubois
*
*/

class ICalFailureModeSvc : virtual public IInterface {
   

public:
   
   
    static const InterfaceID& interfaceID() { return IID_ICalFailureModeSvc; }
    /// get the list of enabled failure mode conditions
    virtual int getFailureConditions()=0;
    /// look for crystal in list of dead towers
    virtual bool matchChannel(idents::CalXtalId id, idents::CalXtalId::XtalFace face)=0;
protected:
  
    /// look for crystal in list of dead towers
    virtual bool matchTower(idents::CalXtalId id)=0;
   
    /// look for crystal in list of dead layers
    virtual bool matchTowerAfee(idents::CalXtalId id, idents::CalXtalId::XtalFace face)=0;
    /// process the input list of towers
    virtual void processTowerList()=0;
    /// process the input list of tower/layer pairs
    virtual void processTowerAfeeList()=0;
};

#endif // ICalFailureModeSvc_H


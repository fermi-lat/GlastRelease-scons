
#ifndef IPropagatorSvc_h
#define IPropagatorSvc_h

#include "GaudiKernel/IInterface.h"

class IKalmanParticle;

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IPropagatorSvc("IPropagatorSvc", 1 , 0); 

/** @class IPropagatorSvc
* @brief Abstract interface for GLAST reconstruction service, giving access to  
* to a IKalmanParticle object.
*
* @author Tracy Usher
* $Header$
*/

class   IPropagatorSvc : virtual public IInterface {
public:
 
    // This returns a pointer to the singleton propagator object
    virtual IKalmanParticle* getPropagator() const = 0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IPropagatorSvc; }
};



#endif  // _H_IPropagatorSvc_

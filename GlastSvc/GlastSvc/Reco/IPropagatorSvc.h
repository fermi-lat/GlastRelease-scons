
#ifndef _H_IPropagatorSvc_
#define _H_IPropagatorSvc_

// includes
//#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/Service.h"
#include "geometry/Ray.h"

#include <memory> // for auto_ptr

//forward declarations
//class IParticlePropagator;
class IKalmanParticle;

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IPropagatorSvc(908, 1 , 0); 

//!  \author: T Burnett
/** Abstract interface for GLAST reconstruction service, giving access to  
   to a IKalmanParticle object.

*/

//class IPropagatorSvc : virtual public IInterface
class IPropagatorSvc 
{
public:
    
    /// return a IParticlePropagator object (auto_ptr guarantees deletion)
    //virtual std::auto_ptr<IParticlePropagator> Propagator(Point p0, Vector v0, double step)=0; 
    //virtual std::auto_ptr<IKalmanParticle> Propagator(Point p0, Vector v0, double step)=0; 

    virtual IKalmanParticle* getPropagator() = 0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IPropagatorSvc; }


};



#endif 

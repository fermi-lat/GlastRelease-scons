// $Id$
// 
//!  \author: Toby Burnett tburnett@u.washington.edu
//

#ifndef _H_IRecoSvc_
#define _H_IRecoSvc_

// includes
#include "GaudiKernel/IInterface.h"
#include "geometry/Point.h"

#include <memory> // for auto_ptr

//forward declarations
class IKalmanParticle;

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IRecoSvc(905, 1 , 0); 

//!  \author: T Burnett
/** Abstract interface for GLAST reconstruction service, giving access to  
   to a IKalmanParticle object.

*/

class   IRecoSvc : virtual public IInterface {
public:
 

    // This returns a pointer to the singleton propagator object
    virtual IKalmanParticle* getPropagator() = 0;


    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IRecoSvc; }

};



#endif  // _H_IGlastDetSvc_

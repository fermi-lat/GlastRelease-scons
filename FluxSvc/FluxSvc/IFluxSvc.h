// $Header$
// 
//!  \author: T. Burnett
//

#ifndef _H_IFluxSvc
#define _H_IFluxSvc

// includes
#include "GaudiKernel/IInterface.h"
#include <string>
#include <list>

// forward declarations
class IFlux;

/** Abstract interface for the flux service

  */
class  IFluxSvc : virtual public IInterface {
public:

    /// just set an IFlux object by name
    virtual  StatusCode source(std::string name, IFlux*&)=0;
    
    /// return a list of legal names
    virtual std::list<std::string> fluxNames()const=0;
};

#endif  // _H_IFluxSvc

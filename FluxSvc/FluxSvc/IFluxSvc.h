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
#include "flux/SpectrumFactoryTable.h"

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IFluxSvc(910, 1 , 0); 

// forward declarations
class IFlux;
class HepRandomEngine;

/** Abstract interface for the flux service

*/
class  IFluxSvc : virtual public IInterface {
public:
    
    /// just get an IFlux object by name
    virtual  StatusCode source(std::string name, IFlux*&)=0;
    
    /// return a list of legal names
    virtual std::list<std::string> fluxNames()const=0;

    /// add a new source
    virtual void addFactory(std::string name, const ISpectrumFactory* factory )=0;
    
    
    /// access to the local HepRandomEngine, to allow synchronization
    virtual HepRandomEngine* getEngine()=0;

    
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IFluxSvc; }
    

};

#endif  // _H_IFluxSvc

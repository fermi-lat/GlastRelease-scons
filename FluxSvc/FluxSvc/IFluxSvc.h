/** 
* @file IFluxSvc.h
* @brief definition of the interface for IFluxSvc
*
*  $Header$
*/
#ifndef _H_IFluxSvc
#define _H_IFluxSvc
/** 
* \class IFluxSvc
*
* \brief The virtual interface for FluxSvc-type objects.
*
* \author Toby Burnett tburnett@u.washington.edu
* 
* $Header $
*/

// includes
#include "GaudiKernel/IInterface.h"
#include "CLHEP/Vector/Rotation.h"

#include <string>
#include <list>
#include <vector>

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IFluxSvc("FluxSvc", 2 , 0); 

// forward declarations
class IFlux;
class IParticlePropertySvc;
class ISpectrumFactory;
class HepRandomEngine;

//! Abstract interface for the flux service, FluxSvc.
class  IFluxSvc : virtual public IInterface {
public:
    
    /// just get an IFlux object by name
    virtual  StatusCode source(std::string name, IFlux*&)=0;
    
    /// return a list of legal names
    virtual std::list<std::string> fluxNames()const=0;
    
    /// add a new source
    virtual void addFactory(std::string name, const ISpectrumFactory* factory )=0;
    
    /// return pointer to the random engine that FluxSvc uses
    virtual HepRandomEngine* getRandomEngine()=0;

    /// pass a specific amount of time
    virtual void pass (double t)=0;    
    
    /// create a set of display windows using rootplot.
    virtual void rootDisplay(std::vector<const char*> arguments)=0;
    
    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IFluxSvc; }
    
    ///return the pointer to the current IFlux object
    virtual IFlux* currentFlux()=0;
    
    /// name of the flux
    virtual std::string fluxName()const=0;
    
    /// set the glast tilt angles.
    virtual void setOrientation(std::pair<double,double> ang)=0;

    /// get the angular values of the satellite
    virtual std::pair<double,double> getOrientation()=0;
    
    
    ///this transforms glast-local (cartesian) vectors into galactic (cartesian) vectors
    virtual HepRotation transformGlastToGalactic(double time)const=0;
    
    /// get the current satellite location
    virtual std::pair<double,double> location()=0;

    virtual void setRockType(int rockType)=0;
    ///this should return the source file names, along with the contained sources.
    virtual std::vector<std::pair< std::string ,std::list<std::string> > > sourceOriginList() const=0;
        
};

#endif  // _H_IFluxSvc

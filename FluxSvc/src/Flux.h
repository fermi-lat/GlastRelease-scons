// $Header$

#ifndef FLUXSVC_FLUX_H
#define FLUXSVC_FLUX_H


#include "FluxSvc/IFlux.h"


// forward declarations
class FluxMgr;
class EventSource;
class FluxSource;


//!  The class holding the interface with FluxMgr, EventSource, and FluxSource of the flux package.
//!  Flux is used to get the actual information(energy, name, etc) about the current particle, and to generate
//!  new ones, through this interface.
class Flux : public IFlux {
public:
    /// ctor, select the name
    Flux(std::string name);
    virtual ~Flux();
    
    /// name of the flux
    virtual std::string name()const;
    
    /// full title of the flux
    virtual std::string title()const;
    
    /// generate a new entry trajectory
    virtual void generate();
    
    /// the particle generated 
    virtual std::string particleName()const;
    
    /// the particle property entry for the last particle generated 
    //virtual ParticleProperty* property()const=0;

    /// its kinetic energy
    virtual double energy()const;
    
    /// starting point 
    virtual HepPoint3D launchPoint()const;
    
    /// direction
    virtual HepVector3D launchDir()const;
    
    /// return the time
    virtual double time()const;

    /// pass a specific amount of time
    virtual void pass ( double t);

    /// Get the time as held by GPS
    /*GPStime*/int gpsTime () const;
    
    /// rate ( /mm**2 /s)
    virtual double rate()const;
    
    /// set the static pointer 
    static void mgr(FluxMgr* );
    
    /// set the area of the target
    virtual void setTargetArea( double area);
    
    /// retrieve the area (a static, same for all fluxes)
    double targetArea()const;
    
    /// find which spectrum created the current particle
    virtual std::string findSource()const;
    
    /// return a unique number correcponding to that spectrum
    virtual int numSource()const;
    
   // virtual void addFactory( const IFactory* factory );
    
    virtual void addFactory(std::string name, const ISpectrumFactory* factory );/* {
                                                                                insert(std::make_pair<std::string, const ISpectrumFactory*>(name, factory));
}*/
    
    
private:
    
    EventSource* m_event;  
    double m_time;  // elapsed time: here for now.
    FluxSource* m_flux; // actual FluxSource used 
    
    static FluxMgr* s_mgr;
    
};

#endif
// $Header$

#ifndef FLUXSVC_FLUX_H
#define FLUXSVC_FLUX_H


#include "FluxSvc/IFlux.h"


// forward declarations
class FluxMgr;
class EventSource;
class FluxSource;

class Flux : public IFlux {
public:
        // ctor, select the name
    Flux(std::string name);
    virtual ~Flux();

    // name of the flux
    virtual std::string name()const;

    /// full title of the flux
    virtual std::string title()const;

    // generate a new entry trajectory
    virtual void generate();

    // the particle generated 
    virtual std::string particleName()const;

    // its kinetic energy
    virtual double energy()const;

    // starting point 
    virtual HepPoint3D launchPoint()const;

    // direction
    virtual HepVector3D launchDir()const;

    // rate ( /mm**2 /s)
    virtual double rate()const;

    // set the static pointer 
    static void mgr(FluxMgr* );
 
    /// set the area of the target
    virtual void setTargetArea( double area);

    /// retrieve the area (a static, same for all fluxes)
    double targetArea()const;

    /// find which spectrum created the current particle
    virtual std::string findSource()const;

    /// return a unique number correcponding to that spectrum
    virtual int numSource()const;

private:

    EventSource* m_event;  
    FluxSource* m_flux; // actual FluxSource used 

    static FluxMgr* s_mgr;

};

#endif
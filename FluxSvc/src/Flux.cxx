// $Header$

// Original author: T. Burnett


#include "Flux.h"

#include "flux/FluxSource.h"
#include "flux/EventSource.h"
#include "flux/FluxMgr.h"

Flux::Flux(std::string name) {
    m_event = s_mgr->source(name);
}
Flux::~Flux() 
{
    delete m_flux;
}

FluxMgr* Flux::s_mgr=0;

void Flux::mgr(FluxMgr* m){ s_mgr=m;}
// name of the flux
std::string Flux::name() const
{
    return m_flux->name();
}

/// full title of the flux
std::string Flux::title()const 
{
    return m_event->fullTitle();
}

// generate a new entry trajectory, set FluxSource
void Flux::generate()
{
    m_flux = m_event->event();
}

// the particle generated 
std::string Flux::particleName()const{
    return m_flux->spectrum()->particleName();
    //TODO: fix for composite spectra
}

// its kinetic energy
double Flux::energy()const
{
    return m_flux->energy();
}

// starting point 
HepPoint3D Flux::launchPoint()const
{
    return m_flux->launchPoint();
}

// direction
HepVector3D Flux::launchDir()const
{
    return m_flux->launchDir();
}

// rate ( /mm**2 /s)
double Flux::rate()const
{
    return m_event->rate();
}

    /// set the area of the target
void Flux::setTargetArea( double area)
{
   m_event->totalArea(area);
}

double Flux::targetArea()const
{
    return m_event->totalArea();
}


/// find which spectrum created the current particle
std::string Flux::findSource()const
{
	return m_event->findSource();
}

/// return a unique number correcponding to that spectrum
int Flux::numSource()const
{
    return m_event->numSource();

}
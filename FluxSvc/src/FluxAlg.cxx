// $Header$

// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"
#include "GaudiKernel/SmartRefVector.h"

// Event for creating the McEvent stuff
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/EventModel.h"

//flux
#include "FluxSvc.h"
#include "FluxSvc/IFlux.h"
#include "Spectrum.h"
#include "SpectrumFactory.h"

#include "CLHEP/Vector/LorentzVector.h"

#include <cassert>
#include <vector>

#include "FluxAlg.h"
//------------------------------------------------------------------------


static const AlgFactory<FluxAlg>  Factory;
const IAlgFactory& FluxAlgFactory = Factory;

//------------------------------------------------------------------------
//! ctor
FluxAlg::FluxAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
{
    // declare properties with setProperties calls
    declareProperty("source_name",  m_source_name="default");
    
}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode FluxAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    if ( service("FluxSvc", m_fluxSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the FluxSvc!" << endreq;
        return StatusCode::FAILURE;
    }
    
    log << MSG::INFO << "loading source..." << endreq;
    
    sc =  m_fluxSvc->source(m_source_name, m_flux);
    if( sc.isFailure()) {
        log << MSG::ERROR << "Could not find flux " << m_source_name << endreq;
        return sc;
    }
    log << MSG::INFO << "Source: "<< m_flux->title() << endreq;
    
    log << MSG::INFO << "Source title: " << m_flux->title() << endreq;
    log << MSG::INFO << "        area: " << m_flux->targetArea() << endreq;
    log << MSG::INFO << "        rate: " << m_flux->rate() << endreq;
    
    if ( service("ParticlePropertySvc", m_partSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the ParticlePropertySvc!" << endreq;
        return StatusCode::FAILURE;
    }
    return sc;
}


//------------------------------------------------------------------------
//! process an event
StatusCode FluxAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    //
    // Purpose: have the flux service create parameters of an incoming particle 
    // if nothing has changed, then use the existing m_flux,
    // but if the "current" IFlux is not the same as the one we have now,
    // then change our m_flux pointer to be the new one.
    // Output:  a staturCode to ensure the function executed properly.
    
    if(m_fluxSvc->currentFlux() == m_flux){
        m_flux->generate();
    }else{
        m_flux = m_fluxSvc->currentFlux();
        m_flux->generate();
    }
    
    HepPoint3D p = m_flux->launchPoint();
    HepPoint3D d = m_flux->launchDir();
    double ke = m_flux->energy(); // kinetic energy in MeV
    std::string particleName = m_flux->particleName();
    //if it's a "timeTick, then ExposureAlg should take care of it, and no othe algorithms should care about it.
    if(particleName == "TimeTick"){
        log << MSG::DEBUG << particleName << " particle found, will exit particle creation/reconstruction loops" << endreq;
        particleName = "gamma"; 
        setFilterPassed( false );//no need to return - we want the time record on the TDS.  the ExposureAlg will handle the rest.
    }
    
    
    //here's where we get the particleID and mass for later.
    if( particleName=="p") particleName="proton";
    ParticleProperty* prop = m_partSvc->find(particleName);
    
    if( prop==0) {
        log << MSG::ERROR << "Particle name " << particleName << " not found by particle properties" << endreq;
        return StatusCode::FAILURE;
    }
    
    int partID = prop->jetsetID(); // same as stdhep id
    
    log << MSG::DEBUG << particleName
        << "(" << m_flux->energy()
        << " MeV), Launch: " 
        << "(" << p.x() <<", "<< p.y() <<", "<<p.z()<<")" 
        << "mm, Dir " 
        << "(" << d.x() <<", "<< d.y() <<", "<<d.z()<<")" 
        << endreq;
    
    
    // Here the TDS is prepared to receive hits vectors
    // Check for the MC branch - it will be created if it is not available
    
    DataObject *mc = new Event::McParticleCol;
    //eventSvc()->retrieveObject("/Event/MC", mc);
    sc=eventSvc()->registerObject(EventModel::MC::Event , mc);
    if(sc.isFailure()) {
        log << MSG::WARNING << EventModel::MC::Event  <<" could not be registered on data store" << endreq;
        //return sc;
    }
    
    
    Event::McParticleCol* pcol = new Event::McParticleCol;
    sc = eventSvc()->registerObject(EventModel::MC::McParticleCol, pcol);
    if( sc.isFailure()) {
        
        log << MSG::ERROR << "Could not Register "<< EventModel::MC::McParticleCol << endreq;
        
        return sc;
    }
    Event::McParticle * parent= new Event::McParticle;
    pcol->push_back(parent);
    
    double mass = prop->mass() , 
        energy = (ke+mass),
        momentum=sqrt(energy*energy - mass*mass); 
    HepLorentzVector pin(d*momentum,energy);
    
    // This parent particle decay at the start in the first particle, 
    // so initial momentum and final one are the same
    parent->initialize(parent, partID, 
        Event::McParticle::PRIMARY,
        pin,p);
    parent->finalize(pin, p);
    
    SmartDataPtr<Event::EventHeader> header(eventSvc(), EventModel::EventHeader);
    if(header) {
        Event::EventHeader& h = header;
        h.setTime(m_flux->time());
        
    } else { 
        log << MSG::WARNING << " could not find the event header" << endreq;
        return StatusCode::SUCCESS;//FAILURE;
    }
    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode FluxAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    return sc;
}


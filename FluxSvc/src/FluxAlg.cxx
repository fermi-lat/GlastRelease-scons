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

// GlastEvent for creating the McEvent stuff
#include "GlastEvent/TopLevel/Event.h"
#include "GlastEvent/TopLevel/MCEvent.h"
#include "GlastEvent/MonteCarlo/McParticle.h"
#include "GlastEvent/TopLevel/EventModel.h"

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

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( FluxAlg );

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

namespace {
    double massOf(std::string name) {
        static double GeV=1.0, MeV = GeV/1000, keV = MeV/1000;
        if( name == "gamma" ) return 0;
        if( name == "p")      return 938.272 * MeV;
        if( name == "e-")     return 511.0 * keV;
        return -1;
    }
}
//------------------------------------------------------------------------
//! process an event
StatusCode FluxAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    //
    // have the flux service create parameters of an incoming particle 
    //
    //using namespace mc;
    
    m_flux->generate();
    
    HepPoint3D p = m_flux->launchPoint();
    HepPoint3D d = m_flux->launchDir();
    double ke = m_flux->energy(); // kinetic energy
    std::string particleName = m_flux->particleName();
    

    //here's where we get the particleID for later.
    ParticleProperty* prop = m_partSvc->find(particleName);

    int partID = prop->type();
      
    log << MSG::DEBUG << particleName
        << "(" << m_flux->energy()
        << " GeV), Launch: " 
        << "(" << p.x() <<", "<< p.y() <<", "<<p.z()<<")" 
        << " Dir " 
        << "(" << d.x() <<", "<< d.y() <<", "<<d.z()<<")" 
        << endreq;
    
    if( m_event%100 ==0 ) log << MSG::WARNING << "Event % " << m_event << endreq;

    mc::McParticleCol* pcol = new mc::McParticleCol;

    // Here the TDS is prepared to receive hits vectors
    // Check for the MC branch - it will be created if it is not available
    DataObject *mc;
    eventSvc()->retrieveObject("/Event/MC", mc);

    //log << MSG::DEBUG << "TDS ready" << endreq;


    IDataProviderSvc* temp = eventSvc();
    //log << MSG::DEBUG << "FluxAlg temp =" << temp << endreq;   
    
    StatusCode sc2 = temp->registerObject("/Event/McParticleCol", pcol);
    if( sc2.isFailure()) {
        log << MSG::ERROR << "Could not Register pcol" << endreq;
        return sc2;
    }
    mc::McParticle * parent= new mc::McParticle;
    pcol->push_back(parent);

    HepLorentzVector pin(d*10.,ke);  //the 10 is so that d is in mm, instead of cm.

    // This parent particle decay at the start in the first particle, 
    // so initial momentum and final one are the same
    parent->initialize(parent, partID, //pdef->GetPDGEncoding(), 
        mc::McParticle::PRIMARY,
        pin);
    parent->finalize(pin, p);
#if 0
    SmartDataPtr<MCEvent> mcEvent(eventSvc(),"/Event/MC");
    if( 0==mcEvent) { 
        log << MSG::ERROR << "could not find the McEvent header" << endreq;
        return StatusCode::FAILURE;
    }
#endif
   // mcEvent->setSourceId(m_flux->numSource());
    
    double mass = massOf(particleName), 
        energy = mass+ke, 
        momentum= sqrt(energy*energy+mass*mass);
    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode FluxAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    return sc;
}

/*
mc::McParticleCol* pcol = new mc::McParticleCol;
    eventSvc()->registerObject("/Event/MC/McParticleCol", pcol);
    mc::McParticle * parent= new mc::McParticle;
    pcol->push_back(parent);

    // This parent particle decay at the start in the first particle, 
    // so initial momentum and final one are the same
    parent->initialize(parent, pdef->GetPDGEncoding(), 
        mc::McParticle::PRIMARY,
        pin);
    parent->finalize(pin, p);
*/
/*// Here the TDS is prepared to receive hits vectors
    // Check for the MC branch - it will be created if it is not available
    DataObject *mc;
    eventSvc()->retrieveObject("/Event/MC", mc);

    log << MSG::DEBUG << "TDS ready" << endreq;*/
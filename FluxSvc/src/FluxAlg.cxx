// $Header$

// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"

// GlastEvent for creating the McEvent stuff
#include "GlastEvent/TopLevel/Event.h"
#include "GlastEvent/TopLevel/MCEvent.h"
#include "GlastEvent/MonteCarlo/McParticle.h"
#include "GlastEvent/MonteCarlo/McVertex.h"

//flux
#include "FluxSvc/FluxSvc.h"
#include "FluxSvc/IFlux.h"
#include "flux/Spectrum.h"
#include "flux/SpectrumFactory.h"

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
    m_flux->generate();
    
    HepPoint3D p = m_flux->launchPoint();
    HepPoint3D d = m_flux->launchDir();
    double ke = m_flux->energy(); // kinetic energy
    std::string particleName = m_flux->particleName();

    ParticleProperty* prop = m_partSvc->find(particleName);

    
    log << MSG::DEBUG << particleName
        << "(" << m_flux->energy()
        << " GeV), Launch: " 
        << "(" << p.x() <<", "<< p.y() <<", "<<p.z()<<")" 
        << " Dir " 
        << "(" << d.x() <<", "<< d.y() <<", "<<d.z()<<")" 
        << endreq;
    
    if( m_event%100 ==0 ) log << MSG::WARNING << "Event % " << m_event << endreq;
    
    SmartDataPtr<MCEvent> mcEvent(eventSvc(),"/Event/MC");
    if( 0==mcEvent) { 
        log << MSG::ERROR << "could not find the McEvent header" << endreq;
        return StatusCode::FAILURE;
    }
    mcEvent->setSourceId(m_flux->numSource());
    
    // set the time (in micro seconds?)
    //event->setTime(static_cast<long>( m_flux->time()*1e6+0.5));
    
    
    // first expect to have converter create this (but empty)
    m_vlist = SmartDataPtr<McVertexCol>(eventSvc(), "/Event/MC/McVertexCol");
    if(m_vlist==0) return StatusCode::FAILURE;
    
    // the converter will make this (also empty)
    DataObject* plist;
    if( (eventSvc()->retrieveObject( "/Event/MC/McParticleCol", plist)).isFailure() ){ 
        return StatusCode::FAILURE;
    }
    try{
        m_plist = dynamic_cast<McParticleCol*>(plist);
    }catch(...) { 
        return StatusCode::FAILURE;
    }
    
    // create the starting node in the Vertex/Particle tree
    McParticle* p1 = new McParticle;
    McVertex*   v1 = new McVertex;
    m_vlist->add(v1);
    m_plist->add(p1);
    
    // a McVertex is really a track segment
    v1->setInitialPosition(p);
    v1->setMcParticle(p1); // associated particle will be set up below
    
    double mass = massOf(particleName), 
        energy = mass+ke, 
        momentum= sqrt(energy*energy+mass*mass);
    v1->setInitialFourMomentum( HepLorentzVector( momentum*d.unit(),energy ) );

    
    // choose among primaryOrigin, daughterOrigin, decayProduct, showerContents, showerBacksplash
    v1->setVertexType(McVertex::primaryOrigin); 
    
    v1->setMotherMcParticle( 0);

#if 0
    p1->setParticleID(p->idCode()); //TODO: is this right?
    p1->setParticleProperty(p->idCode());//TODO: is this right?
    p1->setPrimaryParticleFlag(true);
    p1->setMcVertex(v1);

    
    p1->setParticleID(p->idCode()); //TODO: is this right?
    p1->setParticleProperty(p->idCode());//TODO: is this right?
    p1->setPrimaryParticleFlag(p->status()==MCParticle::PRIMARY);
    p1->setMcVertex(v1);
    if( mother !=0) {
        mother->addDaughterMcParticle(p1);
    } else  m_root = v1;  // make root available for display, etc.
    
    
    HepLorentzVector final;
    for( int i=0; i< p->numChildren(); ++i){
        addParticle(v1, p->child(i));
        final +=  *(p->child(i));
    }
    v1->setFinalFourMomentum(final);
    
#endif
    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode FluxAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    return sc;
}


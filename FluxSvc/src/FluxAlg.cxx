/** @file FluxAlg.cxx
@brief declaration and definition of the class FluxAlg

$Header$

*/

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

#include "EventSource.h"

#include "CLHEP/Vector/LorentzVector.h"

#include <cassert>
#include <vector>



// Include files
// Gaudi system includes
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"

class IFlux;
class IFluxSvc;
class IparticlePropertySvc;


/** 
* \class FluxAlg
*
* \brief This is an Algorithm designed to get particle information 
* from FluxSvc and put it onto the TDS for later retrieval
* \author Toby Burnett
* 
* $Header$
*/

class FluxAlg : public Algorithm {
public:
    FluxAlg(const std::string& name, ISvcLocator* pSvcLocator);
    double currentRate(){return m_currentRate;}
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
    
private: 
    double m_currentRate;
    StringProperty m_source_name;

    IFluxSvc*   m_fluxSvc;
    IFlux *     m_flux;
    
    
    UnsignedIntegerProperty m_run;      // run number
    unsigned int m_sequence;  // sequence number
    
    
    IDataProviderSvc* m_eds;
    
    IParticlePropertySvc * m_partSvc;

    // the target area
    DoubleProperty m_area;
    IntegerProperty m_pointing_mode;
    DoubleProperty m_rocking_angle; // x-axis
    DoubleProperty m_rocking_angle_z; // z-axis

    std::map<int, int> m_counts; //! for measuring the total number generated per code.
    
};
//------------------------------------------------------------------------


static const AlgFactory<FluxAlg>  Factory;
const IAlgFactory& FluxAlgFactory = Factory;

//------------------------------------------------------------------------
//! ctor
FluxAlg::FluxAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator) , m_sequence(0)
{
    // declare properties with setProperties calls
    declareProperty("source_name",  m_source_name="default");
    declareProperty("MCrun",        m_run=100);
    declareProperty("area",        m_area=6.0); // target area in m^2
    declareProperty("pointing_mode", m_pointing_mode=0);
    declareProperty("rocking_angle", m_rocking_angle=0); // in degrees
    declareProperty("rocking_angle_z", m_rocking_angle_z=0); // in degrees
    
}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode FluxAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    // set target area for random point generation
    EventSource::totalArea(m_area);
    
    // set pointing mode and associated rocking angle 
    if ( service("FluxSvc", m_fluxSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the FluxSvc!" << endreq;
        return StatusCode::FAILURE;
    }
    //this line sets the explicit rocking angles to be used IF the 
    //rocking type is "explicit." Note that it uses the m_rocking_angle also.
    m_fluxSvc->setExplicitRockingAngles(m_rocking_angle*M_PI/180,m_rocking_angle_z*M_PI/180);
    //m_fluxSvc->setRockingAngle(m_rocking_angle);

    //then this line sets the rocking type, as well as the rocking angle.
    m_fluxSvc->setRockType(m_pointing_mode,m_rocking_angle);

    log << MSG::INFO << "loading source..." << endreq;
    
    sc =  m_fluxSvc->source(m_source_name, m_flux);
    if( sc.isFailure()) {
        log << MSG::ERROR << "Could not find flux " << m_source_name << endreq;
        return sc;
    }
   
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
    HepVector3D d = m_flux->launchDir();
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
        << " mm, Dir " 
        << "(" << d.x() <<", "<< d.y() <<", "<<d.z()<<")" 
        << endreq;
    
    
    // Here the TDS is prepared to receive hits vectors
    // Check for the MC branch - it will be created if it is not available
    Event::MCEvent* mch = 0;

    SmartDataPtr<Event::MCEvent> mcheader(eventSvc(), EventModel::MC::Event);
    if (mcheader == 0) {
        sc=eventSvc()->registerObject(EventModel::MC::Event , mch= new Event::MCEvent);
        mch->initialize(0,0,m_sequence);
        if(sc.isFailure()) {
            log << MSG::WARNING << EventModel::MC::Event  <<" could not be registered on data store" << endreq;
            return sc;
        }
        
    }else mch = mcheader;


    mch->initialize(mch->getRunNumber(), m_flux->numSource(), mch->getSequence());


    
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
    
    // get the event header to set the time
    Event::EventHeader* h = 0; 

    SmartDataPtr<Event::EventHeader> header(eventSvc(), EventModel::EventHeader);
    if(0==header) {
        // not already there: try to register instead
        sc = eventSvc()->registerObject(EventModel::EventHeader, h=new Event::EventHeader);
        if( sc.isFailure()) {
            log << MSG::WARNING << " could not find or register the event header" << endreq;
        }
    } else{  h = header;
    }

    h->setTime(m_flux->time());

    int numEvents = ++m_sequence;
    m_currentRate=numEvents/(m_flux->time());
    m_counts[m_flux->numSource()]++; // update count fo r
    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode FluxAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "Computed Rate: "<< currentRate() << " Hz"
        << "\n\t\t\tSource ID    counts";
    for(std::map<int,int>::const_iterator im=m_counts.begin(); im !=m_counts.end(); ++im) {
        log << "\n\t\t\t\t" << im->first << "\t" << im->second;
    }
    log  << endreq;
    
    return sc;
}


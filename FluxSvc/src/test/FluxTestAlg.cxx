// $Header$

// Include files
#include "FluxSvc/IFluxSvc.h"
#include "FluxSvc/IFlux.h"


// Event for creating the McEvent stuff
//#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/EventModel.h"

// Gaudi system includes
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include <list>
#include <string>
#include <vector>
#include "GaudiKernel/ParticleProperty.h"

//#include "FluxAlg.h"
/*! \class FluxTestAlg
\brief 
In addition to the normal Gaudi JobOptions requirements, there are:
FluxSvc.source_lib, which should contain the relevant xml files to be used. and
FluxTestAlg.source_name, which holds the name of the desired spectrum.
*/
//class ParticleProperty;

class FluxTestAlg : public Algorithm {
    
public:
    //! Constructor of this form must be provided
    FluxTestAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
    
private:
    IFlux* m_flux;
    std::string m_source_name;
    IParticlePropertySvc * m_partSvc;
};


static const AlgFactory<FluxTestAlg>  Factory;
const IAlgFactory& FluxTestAlgFactory = Factory;
/*
void WARNING (const char * text ){  std::cerr << "WARNING: " << text << '\n';}
void FATAL(const char* s){std::cerr << "\nERROR: "<< s;}
*/

//------------------------------------------------------------------------------
//
FluxTestAlg::FluxTestAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator){
    
    declareProperty("source_name", m_source_name="default");
}


//------------------------------------------------------------------------------
/*! */
StatusCode FluxTestAlg::initialize() {
    
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initializing..." << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    // get the pointer to the flux Service 
    IFluxSvc* fsvc;
    
    // get the service
    StatusCode sc = service("FluxSvc", fsvc);
    
    if( sc.isFailure()) {
        log << MSG::ERROR << "Could not find FluxSvc" << endreq;
        return sc;
    }
    
    
    log << MSG::INFO << "loading source..." << endreq;
    
    
    sc =  fsvc->source(m_source_name, m_flux);
    if( sc.isFailure()) {
        log << MSG::ERROR << "Could not find flux " << m_source_name << endreq;
        return sc;
    }
    
    // then do the output here.
    log << MSG::INFO << "start of other loops" << endreq;
    log << MSG::INFO << "Source title: " << m_flux->title() << endreq;
    log << MSG::INFO << "       area: " << m_flux->targetArea() << endreq;
    log << MSG::INFO << "       rate: " << m_flux->rate() << endreq;
    
    
    
    return sc;
}


//------------------------------------------------------------------------------
StatusCode FluxTestAlg::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );    
    
    //Event::McParticleCol*  pcol2= SmartDataPtr<Event::McParticleCol>(eventSvc(), "/Event/MC/McParticleCol");
    
    Event::McParticleCol* pcol = new Event::McParticleCol;
    eventSvc()->retrieveObject("/Event/MC/McParticleCol",(DataObject *&)pcol);
    
    HepVector3D p,d;
    double energy;
    std::string partName;
    //only make a new source if one does not already exist.
    if(pcol==0){
        m_flux->generate();
        p = m_flux->launchPoint()*10.;
        d = m_flux->launchDir();
        energy = m_flux->energy();
        partName = m_flux->particleName();
    }else{
        Event::McParticleCol::iterator elem = (*pcol).begin();
        d = (*elem)->initialFourMomentum().v()*10;
        p = (*elem)->finalPosition();
        
        energy = (*elem)->initialFourMomentum().e();
        /*StdHepId*/int pID = (*elem)->particleProperty();
        if ( service("ParticlePropertySvc", m_partSvc).isFailure() ){
            log << MSG::ERROR << "Couldn't find the ParticlePropertySvc!" << endreq;
            return StatusCode::FAILURE;
        }
        partName = m_partSvc->findByStdHepID(pID)->particle();
    }   
    
    log << MSG::INFO << partName
        << "(" << energy
        << " GeV), Launch: " 
        << "(" << p.x() <<", "<< p.y() <<", "<<p.z()<<")" 
        << " Dir " 
        << "(" << d.x() <<", "<< d.y() <<", "<<d.z()<<")"
         << ",  Elapsed Time = " << m_flux->time()
        << endreq;
    
    //m_flux->pass(10.);
    return sc;
}


//------------------------------------------------------------------------------
StatusCode FluxTestAlg::finalize() {
    
    return StatusCode::SUCCESS;
}








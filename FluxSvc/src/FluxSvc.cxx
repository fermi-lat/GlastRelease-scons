// $Header$
// 
//  Original author: Toby Burnett tburnett@u.washington.edu
//

#include "FluxSvc/FluxSvc.h"

#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/IRndmGenSvc.h"

#include "GaudiKernel/IParticlePropertySvc.h"

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RanluxEngine.h"


#include "Flux.h"

#include "flux/FluxMgr.h"
#include <algorithm>

// declare the service factories for the FluxSvc
static SvcFactory<FluxSvc> a_factory;
const ISvcFactory& FluxSvcFactory = a_factory;

//void FATAL(const char* text){std::cerr << text << std::endl;}
//void WARNING(const char* text){std::cerr << text << std::endl;}

static std::string default_source_library("$(FLUXROOT)/xml/source_library.xml");

// ------------------------------------------------
// Implementation of the FluxSvc class
// ------------------------------------------------
/// Standard Constructor
FluxSvc::FluxSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc)
{
    
    declareProperty("source_lib" , m_source_lib);   
  
    HepRandom::setTheEngine(new RanluxEngine);
}

std::list<std::string> FluxSvc::fluxNames()const{
    return m_fluxMgr->sourceList();
}

StatusCode FluxSvc::source(std::string name, IFlux*& flux) {
    
    std::list<std::string> source_list( fluxNames() );
    std::list<std::string> source_list2( SpectrumFactoryTable::instance()->spectrumList() );
    
    if( std::find(source_list.begin(), source_list.end(), name) == source_list.end() 
        &&(std::find(source_list2.begin(), source_list2.end(), name) == source_list2.end()))
        //flux =  new Flux(name);
        return StatusCode::FAILURE;
    
    flux =  new Flux(name);
    
    return StatusCode::SUCCESS;
}

/// Standard Destructor
FluxSvc::~FluxSvc()  
{
}


// initialize
StatusCode FluxSvc::initialize () 
{
    //std::vector<std::string> fileList;
    
    StatusCode  status =  Service::initialize ();
    
    // bind all of the properties for this service
    setProperties ();
    
    // open the message log
    MsgStream log( msgSvc(), name() );

    //If source library was not set, put in default
    if( m_source_lib.empty() ){
        m_source_lib.push_back(default_source_library);
        log << MSG::INFO << "Set source library list to " << default_source_library << endreq;
    }

     
    // create a FluxMgr object which will then be available.
    m_fluxMgr = new FluxMgr(m_source_lib);
    
    Flux::mgr(m_fluxMgr); // tell our Flux object
    
    // check that it was made properly
    if( m_fluxMgr->sourceList().empty()) {
        log << MSG::ERROR  << "Did not initialize properly: no sources detected" << endreq;
        status = StatusCode::FAILURE;
    }
        
    if ( service("ParticlePropertySvc", m_partSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the ParticlePropertySvc!" << endreq;
        return StatusCode::FAILURE;
    }

    
    return status;
}

// finalize
StatusCode FluxSvc::finalize ()
{
    StatusCode  status = StatusCode::SUCCESS;
    
    delete m_fluxMgr;
    return status;
}

/// Query interface
StatusCode FluxSvc::queryInterface(const IID& riid, void** ppvInterface)  {
    if ( IID_IFluxSvc.versionMatch(riid) )  {
        *ppvInterface = (IFluxSvc*)this;
    }
    else  {
        return Service::queryInterface(riid, ppvInterface);
    }
    addRef();
    return SUCCESS;
}


/// add a new source
void FluxSvc::addFactory(std::string name, const ISpectrumFactory* factory ){
    m_fluxMgr->addFactory(name, factory);
}

/// access to the local HepRandomEngine, to allow synchronization
HepRandomEngine* FluxSvc::getEngine()
{
    return HepRandom::getTheEngine();
}
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

// nasty include needed so that we can get at the HepRandomEngine
// in the initialize method. (Ian Gable) 
// The actual is GaudiSvc
#if 0
#include "src/RndmGenSvc/HepRndmBaseEngine.h"
#endif
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RanluxEngine.h"


#include "Flux.h"

#include "flux/FluxMgr.h"
#include <algorithm>

// declare the service factories for the FluxSvc
static SvcFactory<FluxSvc> a_factory;
const ISvcFactory& FluxSvcFactory = a_factory;

void FATAL(const char* text){std::cerr << text << std::endl;}
void WARNING(const char* text){std::cerr << text << std::endl;}


// ------------------------------------------------
// Implementation of the FluxSvc class
// ------------------------------------------------
/// Standard Constructor
FluxSvc::FluxSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc)
{
    // declare the properties and set defaults
    
    //declareProperty ("size", m_default_name);
	
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
    StatusCode  status =  Service::initialize ();
        
    // bind all of the properties for this service
    setProperties ();
    
    // open the message log
    MsgStream log( msgSvc(), name() );
 
    // create a FluxMgr object which will then be available.
    m_fluxMgr = new FluxMgr();

    Flux::mgr(m_fluxMgr); // tell our Flux object

    // check that it was made properly
    if( m_fluxMgr->sourceList().empty()) {
        log << MSG::ERROR  << "Did not initialize properly: no sources detected" << endreq;
        status = StatusCode::FAILURE;
    }

#if 0
    // set up the random number service to use Gaudi's default (in progress )
    IRndmGenSvc* randSvc = 0;
    StatusCode sc = service( "RndmGenSvc", randSvc, true );
    if( sc.isFailure() ) { log << MSG::ERROR << "Failed to get random service " << endreq; return sc; }

    IRndmEngine* engine = 0;
    if( (engine=randSvc->engine())==0 ) {
        log << MSG::ERROR << "Failed to get the IRndmEngine" << endreq;
    }
   

    // cast to the base type so that we can actually get access to the 
    // HepRandomEngine.
    HepRndm::BaseEngine* baseEngine = dynamic_cast<HepRndm::BaseEngine*>(engine);
    
    if(!baseEngine)
    {
        log << MSG::ERROR << "Failed to cast IRndmEngine to BaseEngine" <<endreq;
        return StatusCode::FAILURE;
    }
    
    
    HepRandomEngine* hepengine = baseEngine->hepEngine();
    
    if(!hepengine)
    {
        log << MSG::ERROR << "No HepRandomEngine available from gaudi!" << endreq;
        return StatusCode::FAILURE;
    }
    
    HepRandom* gen = HepRandom::getTheGenerator(); 
    gen->setTheEngine(hepengine);
   
#endif

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

HepRandomEngine* FluxSvc::getEngine()
{
	return HepRandom::getTheEngine();
}
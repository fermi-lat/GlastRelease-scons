// $Header$
// 
//  Original author: Toby Burnett tburnett@u.washington.edu
//

#include "FluxSvc.h"
#include "FluxSvc/IRegisterSource.h"


#include "./test/flux/rootplot.h"

#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/IObjManager.h"
#include "GaudiKernel/IToolFactory.h"

#include "GaudiKernel/IParticlePropertySvc.h"
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RanluxEngine.h"

#include "Flux.h"

#include "FluxMgr.h"
#include <algorithm>

// declare the service factories for the FluxSvc
static SvcFactory<FluxSvc> a_factory;
const ISvcFactory& FluxSvcFactory = a_factory;


static std::string default_source_library("$(FLUXSVCROOT)/xml/source_library.xml");
static std::string default_dtd_file("$(FLUXSVCROOT)/xml/source.dtd");

// ------------------------------------------------
// Implementation of the FluxSvc class
// ------------------------------------------------
/// Standard Constructor
FluxSvc::FluxSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc), m_currentFlux(0)
{
    
    declareProperty("source_lib" , m_source_lib); 
    declareProperty("dtd_file"   , m_dtd_file=default_dtd_file);
    
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
    m_currentFlux = flux;    
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
    
    // If source library was not set, put in default
    if( m_source_lib.empty() ){
        m_source_lib.push_back(default_source_library);
        log << MSG::INFO << "Set source library list to " << default_source_library << endreq;
    }
    
    
    // create a FluxMgr object which will then be available.
    m_fluxMgr = new FluxMgr(m_source_lib, m_dtd_file);
    
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
    
    //----------------------------------------------------------------
    // most of  the following cribbed from ToolSvc and ObjManager
    
    // look for a factory of an AlgTool that implements the IRegisterSource interface:
    // if found, make one and call the special method 
    
    // Manager of the AlgTool Objects
    IObjManager* objManager=0;             
    
    // locate Object Manager to locate later the tools 
    status = serviceLocator()->service("ApplicationMgr", objManager );
    if( status.isFailure()) {
        log << MSG::ERROR << "Unable to locate ObjectManager Service" << endreq;
        return status;
    }
    
    
    IToolSvc* tsvc  =0;
    status = service( "ToolSvc", tsvc, true );
    if( status.isFailure() ) {
        log << MSG::ERROR << "Unable to locate Tool Service" << endreq;
        return status;
    }
    
    IToolFactory* toolfactory = 0;
    
    // search throught all objects (factories?)
    for(IObjManager::ObjIterator it = objManager->objBegin(); it !=objManager->objEnd(); ++ it){
        
        std::string tooltype= (*it)->ident();
        // is it a tool factory?
        const IFactory* factory = objManager->objFactory( tooltype );
        IFactory* fact = const_cast<IFactory*>(factory);
        status = fact->queryInterface( IID_IToolFactory, (void**)&toolfactory );
        if( status.isSuccess() ) {
            
            // yes: now see if the tool implements the IRegisterSource interface
            IRegisterSource* ireg;
            status = tsvc->retrieveTool(tooltype, ireg);
            if( status.isSuccess() ){
                log << MSG::DEBUG << "Registering sources in " << tooltype << endreq;
                ireg->registerMe(this);
                tsvc->releaseTool(ireg);
            }
            
        }
        
    }
    
    
    return StatusCode::SUCCESS;
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

/// pass a specific amount of time
void FluxSvc::pass ( double t){
    m_fluxMgr->pass(t);
}

void FluxSvc::rootDisplay(std::vector<char*> arguments){
    rootplot abc(arguments);
}

///return the pointer to the current IFlux object
IFlux* FluxSvc::currentFlux(){
    return m_currentFlux;
}

/// name of the flux
std::string FluxSvc::fluxName()const{
    return m_currentFlux->name();
}

void FluxSvc::setOrientation(std::pair<double,double> ang){
    m_fluxMgr->setOrientation(ang);
}

/// get the angular values of the satellite
std::pair<double,double> FluxSvc::getOrientation(){
    return m_fluxMgr->getOrientation();
}

Rotation FluxSvc::transformGlastToGalactic(double time)const{
    return m_fluxMgr->transformGlastToGalactic(time);
}

std::pair<double,double> FluxSvc::location(){
    return m_fluxMgr->location();
}

void WARNING (const char * text ){  std::cerr << "WARNING: " << text << '\n';}
void FATAL(const char* s){std::cerr << "\nERROR: "<< s;}

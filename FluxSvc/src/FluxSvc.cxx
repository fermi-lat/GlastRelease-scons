/** 
* @file FluxSvc.cxx
* @brief definition of the class FluxSvc
*
*  $Header$
*  Original author: Toby Burnett tburnett@u.washington.edu
*/

#include "FluxSvc.h"
#include "FluxSvc/IRegisterSource.h"


#include "FluxSvc/../src/test/flux/rootplot.h"

#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/IObjManager.h"
#include "GaudiKernel/IToolFactory.h"
#include "GaudiKernel/IAlgManager.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/IParticlePropertySvc.h"

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
    declareProperty("EvtMax"     , m_evtMax=0);
    declareProperty("StartTime"   , m_startTime=0);
    declareProperty("EndTime",      m_endTime=0);
    
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
    
    status = serviceLocator()->queryInterface(IID_IAppMgrUI, (void**)&m_appMgrUI);

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

            IAlgTool* itool;
            status = tsvc->retrieveTool(tooltype, itool);
            if( status.isSuccess()) { 
                status =itool->queryInterface( IRegisterSource::interfaceID(), (void**)&itool);
                if( status.isSuccess() ){
                    log << MSG::INFO << "Registering sources in " << tooltype << endreq;
                    dynamic_cast<IRegisterSource*>(itool)->registerMe(this);
                }else{
		  tsvc->releaseTool(itool);
		}

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
    }else if (IID_IRunable.versionMatch(riid) ) {
      *ppvInterface = (IRunable*)this;
    } else  {
        return Service::queryInterface(riid, ppvInterface);
    }

    addRef();
    return SUCCESS;
}


/// add a new source
void FluxSvc::addFactory(std::string name, const ISpectrumFactory* factory ){
    m_fluxMgr->addFactory(name, factory);
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

/// return a string which uniquely identifies the source
std::string FluxSvc::uniqueIDString()const{
    std::strstream t;
    t << m_currentFlux->numSource();
    return m_currentFlux->name() + t.str();
}


void FluxSvc::setOrientation(std::pair<double,double> ang){
    m_fluxMgr->setOrientation(ang);
}

/// get the angular values of the satellite
std::pair<double,double> FluxSvc::getOrientation(){
    return m_fluxMgr->getOrientation();
}

HepRotation FluxSvc::transformGlastToGalactic(double time)const{
    return m_fluxMgr->transformGlastToGalactic(time);
}

std::pair<double,double> FluxSvc::location(){
    return m_fluxMgr->location();
}
/// this sets the rocking mode in GPS.
void FluxSvc::setRockType(int rockType){
   m_fluxMgr->setRockType(rockType);
}

std::vector<std::pair< std::string ,std::list<std::string> > > FluxSvc::sourceOriginList() const{
    return m_fluxMgr->sourceOriginList();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode FluxSvc::run(){
    StatusCode status = StatusCode::FAILURE;
    MsgStream log( msgSvc(), name() );

    if ( 0 == m_appMgrUI )  return status; 

    IProperty* propMgr=0;
    status = serviceLocator()->service("ApplicationMgr", propMgr );
    if( status.isFailure()) {
        log << MSG::ERROR << "Unable to locate PropertyManager Service" << endreq;
        return status;
    }
    
    IntegerProperty evtMax("EvtMax",0);
    status = propMgr->getProperty( &evtMax );
    if (status.isFailure()) return status;
    
    setProperty(evtMax);
    
    // now find the top alg so we can monitor its error count
    //
    IAlgManager* theAlgMgr;
    status = serviceLocator( )->getService( "ApplicationMgr",
        IID_IAlgManager,
        (IInterface*&)theAlgMgr );
    IAlgorithm* theIAlg;
    Algorithm*  theAlgorithm=0;
    IntegerProperty errorProperty("ErrorCount",0);
    
    status = theAlgMgr->getAlgorithm( "Top", theIAlg );
    if ( status.isSuccess( ) ) {
        try{
            theAlgorithm = dynamic_cast<Algorithm*>(theIAlg);
        } catch(...){
            status = StatusCode::FAILURE;
        }
    }
    if ( status.isFailure( ) ) {
        log << MSG::WARNING << "Could not find algorithm 'Top'; will not monitor errors" << endreq;
    }
    
    
    // loop over the events
    IFlux* flux=currentFlux();
    int eventNumber= 0;
    double currentTime=m_startTime;
    
    { bool noend=true;
    log << MSG::INFO << "Runable interface starting event loop as :" ; 
    if( m_evtMax>0)  { log << " MaxEvt = " << m_evtMax; noend=false;  }
    if( m_endTime>0) { log << " EndTime= " << m_endTime; noend=false; }
    log << endreq;
    
    if(noend) { 
        log << MSG::WARNING << "No end condition specified: will not process any events!" << endreq; 
    }
    }
    while( m_evtMax>0 && eventNumber < m_evtMax
        || m_endTime>0 && currentTime< m_endTime ) {
        
        status =  m_appMgrUI->nextEvent(1); // currently, always success
        
        // the single event may have created a failure. Check the ErrorCount propery of the Top alg.
        if( theAlgorithm !=0) theAlgorithm->getProperty(&errorProperty);
        if( status.isFailure() || errorProperty.value() > 0){
            status = StatusCode::FAILURE;
        }
        
        if( status.isFailure()) break;
        currentTime = flux->gpsTime();
        eventNumber ++;
    }
    if( status.isFailure()){
        log << MSG::ERROR << "Terminating FluxSvc loop due to error" << endreq;
        
    }else if( m_endTime>0 && currentTime >= m_endTime ) {
        log << MSG::INFO << "Loop terminated by time " << endreq;
    }else {
        log << MSG::INFO << "Processing loop terminated by event count" << endreq;
    }
    return status;
    
}


// Implementation file for GlastRandomSvc which gets adresses 
// of CLHEP random number engines used in Gleam shared libraries
// and sets seeds for them based on run and particle sequence
// number obtained from the MCHeader
//
// $Header$
//
// Author: Toby Burnett, Karl Young


#include "GlastRandomSvc.h"

// #include <map>
#include <vector>
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataStoreItem.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/IObjManager.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/IToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "GlastSvc/GlastRandomSvc/IRandomAccess.h"
#include "GlastSvc/GlastRandomSvc/RandomAccess.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/TopLevel/EventModel.h"

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RandomEngine.h"


// declare the service factories for the GlastRandomSvc
static SvcFactory<GlastRandomSvc> a_factory;
const ISvcFactory& GlastRandomSvcFactory = a_factory; 

GlastRandomSvc::GlastRandomSvc(const std::string& name,ISvcLocator* svc) : Service(name,svc)
{
    // declare the properties
    declareProperty("RandomEngine",  m_randomEngine="TripleRand");
    declareProperty("RunNumber",      m_RunNumber=10);
    declareProperty("InitialSequenceNumber", m_InitialSequenceNumber=0);
}

GlastRandomSvc::~GlastRandomSvc()  
{
    // this is  repetetive in case finalize was
    // not called.
    finalize();
}


StatusCode  GlastRandomSvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_IGlastRandomSvc == riid) {
        *ppvIF = dynamic_cast<IGlastRandomSvc*> (this);
        return StatusCode::SUCCESS;
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
}


StatusCode GlastRandomSvc::initialize () 
{
    StatusCode  status = StatusCode::SUCCESS;
    
    // Open the message log
    MsgStream log( msgSvc(), name() );
    
    // Call super-class
    Service::initialize ();
    
    // Bind all of the properties for this service
    if ( (status = setProperties()).isFailure() ) {
        log << MSG::ERROR << "Failed to set properties" << endreq;
    }
    m_SequenceNumber = m_InitialSequenceNumber;

    // use the incident service to register begin, end events
    IIncidentSvc* incsvc = 0;
    StatusCode sc = service ("IncidentSvc", incsvc, true);
    
    if( status.isFailure() ) return status;
    
    incsvc->addListener(this, "BeginEvent", 0);
    
    // Look for a factory of an AlgTool that implements the
    // IRandomAccess interface:
    // if found, make one and call the special method
    // Manager of the AlgTool Objects
    IObjManager* objManager=0;
    
    // locate Object Manager to locate later the tools
    status = serviceLocator()->service("ApplicationMgr",
        objManager );
    if( status.isFailure()) {
        log << MSG::ERROR << "Unable to locate ObjectManager Service"
            << endreq;
        return status;
    }
    
    IToolSvc* tsvc  =0;
    status = service( "ToolSvc", tsvc, true );
    if( status.isFailure() ) {
        log << MSG::ERROR << "Unable to locate Tool Service" << endreq;
        return status;
    }
    
    IToolFactory* toolfactory = 0;
    
    // Search throught all objects (factories?)
    for(IObjManager::ObjIterator it = objManager->objBegin(); it
        !=objManager->objEnd(); ++ it){
        log << MSG::DEBUG << "Checking object " << it << endreq;
        std::string tooltype= (*it)->ident();
        log << MSG::DEBUG << "of type " << tooltype << endreq;
        // is it a tool factory?
        const IFactory* factory = objManager->objFactory( tooltype );
        log << MSG::DEBUG << "Is " << factory << " a tool factory ?" << endreq;
        IFactory* fact = const_cast<IFactory*>(factory);
        status = fact->queryInterface( IID_IToolFactory,
            (void**)&toolfactory );
        if( status.isSuccess() ) {
            log << MSG::DEBUG << "Guess so... is it a RandomAccess instance ?" << endreq; 
            IAlgTool* itool;
            status = tsvc->retrieveTool(tooltype, itool);
            if( status.isSuccess()) {
                log << MSG::DEBUG << "Ok retrieved it..." << endreq;
                status =itool->queryInterface(IRandomAccess::interfaceID(), 
                    (void**)&itool);
            }
            if( status.isSuccess() ){
                log << MSG::DEBUG << "Found random access: " << tooltype 
                    << ", setting engine to " << m_randomEngine << endreq;
                // Set the Random engine
                dynamic_cast<IRandomAccess*>(itool)->setTheEngine(m_randomEngine);
                // Get its address
                HepRandomEngine* hr =
                    dynamic_cast<IRandomAccess*>(itool)->getTheEngine();
                log << MSG::DEBUG << "Random engine at: " << (int)hr << endreq;
                // Store its name and address in a map
                m_engineMap[tooltype] = hr;
            }
            tsvc->releaseTool(itool);
        }
    }
    return StatusCode::SUCCESS;
}

StatusCode GlastRandomSvc::finalize () {return StatusCode::SUCCESS;}

const IID&  GlastRandomSvc::type () const {
    return IID_IGlastRandomSvc;
}

// handle "incidents" - i.e. set seeds at beginning of event loop
void GlastRandomSvc::handle(const Incident &inc)
{
    
    // Open the message log
    MsgStream log( msgSvc(), name() );
    
    if( inc.type()=="BeginEvent") {
        IDataProviderSvc* eventSvc;
        StatusCode status;
        if( serviceLocator() ) {
            status = serviceLocator()->service("EventDataSvc", eventSvc, true );
        }
        if(status.isFailure())
        {
            log << MSG::ERROR << "Could not find EventDataSvc" << endreq;
        }
        // See if MCEvent was set up properly
        SmartDataPtr<Event::MCEvent> mcevt(eventSvc, EventModel::MC::Event);
        if (mcevt == 0) {
            log << MSG::ERROR << "Error accessing MCEvent" << endreq;
            return;
        }

        mcevt->initialize(m_RunNumber, -1, m_SequenceNumber++);
        int multiplier = 1; 
        int dummy = 0; // for 2nd argument to setSeed
        std::map< std::string, HepRandomEngine* >::const_iterator dllEngine;
        for (dllEngine = m_engineMap.begin(); dllEngine != m_engineMap.end(); ++dllEngine ) {
            long theSeed = multiplier * 100000 * ((m_RunNumber+1) % 20000) + ((m_SequenceNumber+1) % 100000);
            log << MSG::DEBUG << "Setting seed for " <<  dllEngine->first << " to " <<  theSeed << endreq;
            dllEngine->second->setSeed(theSeed,dummy);
            ++multiplier;
        }
    }
}

// Implementation file for GlastRandomSvc which gets adresses 
// of CLHEP random number engines used in Gleam shared libraries
// and sets seeds for them based on run and particle sequence
// number obtained from the MCHeader
//
// $Header$
//
// Author: Toby Burnett, Karl Young


#include "GlastRandomSvc.h"

#include <iterator>
#include <fstream>
#include <iterator>
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
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RandomEngine.h"


// declare the service factories for the GlastRandomSvc
static SvcFactory<GlastRandomSvc> a_factory;
const ISvcFactory& GlastRandomSvcFactory = a_factory; 

GlastRandomSvc::GlastRandomSvc(const std::string& name,ISvcLocator* svc) : Service(name,svc)
{
  // Purpose and Method: Constructor - Declares and sets default properties
  //                     
  // Inputs: service name and locator 
  //         
  // Outputs: None
  // Dependencies: None
  // Restrictions and Caveats:  None

  // declare the properties
  declareProperty("RandomEngine",  m_randomEngine="TripleRand");
  declareProperty("RunNumber",      m_RunNumber=10);
  declareProperty("InitialSequenceNumber", m_InitialSequenceNumber=0);
  declareProperty("SeedFile", m_seedFile="");
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
    status = service ("IncidentSvc", incsvc, true);
    
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
    
  // read seeds from file
  if(m_seedFile != "") {

    std::ifstream input(m_seedFile.c_str());
    if(!input) {
      log << MSG::ERROR << "Unable to find seed file: " << m_seedFile.c_str()
	  << endreq;
      status = StatusCode::FAILURE;
      return status;
    }

    std::istream_iterator<GlastRandomSeed> begin(input);
    std::istream_iterator<GlastRandomSeed> end;
    std::copy(begin, end, std::back_inserter(m_seeds));

    // set no of events
    IProperty* propMgr=0;
    status = serviceLocator()->service("ApplicationMgr", propMgr );
    if( status.isFailure()) {
        log << MSG::ERROR << 
	  "Unable to locate PropertyManager Service" << endreq;
        return status;
    }
    
    IntegerProperty evtMax("EvtMax", m_seeds.size());
    status = propMgr->setProperty(evtMax);
    if (status.isFailure()) return status;
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
        std::string tooltype= (*it)->ident();
        // is it a tool factory?
        const IFactory* factory = objManager->objFactory( tooltype );
        IFactory* fact = const_cast<IFactory*>(factory);
        status = fact->queryInterface( IID_IToolFactory,
            (void**)&toolfactory );
        if( status.isSuccess() ) {
            IAlgTool* itool;
            status = tsvc->retrieveTool(tooltype, itool);
            if( status.isSuccess()) {
                status =itool->queryInterface(IRandomAccess::interfaceID(), 
                    (void**)&itool);
            }
            if( status.isSuccess() ){
	      // note: the following should be DEBUG rather than INFO
	      //       but for some reason DEBUG won't print here
	      // log << MSG::INFO << "Found RandomAccess tool: " << itool << endreq;
                // Set the Random engine
                dynamic_cast<IRandomAccess*>(itool)->setTheEngine(m_randomEngine);
                // Get its address
                HepRandomEngine* hr =
                    dynamic_cast<IRandomAccess*>(itool)->getTheEngine();
		// note: the following should be DEBUG rather than INFO
		//       but for some reason DEBUG won't print here
		// log << MSG::INFO << "CLHEP Engine for " << itool << " at " << hr << endreq;
		// Store its name and address in a map
                m_engineMap[tooltype] = hr;
            }else{
            tsvc->releaseTool(itool);
	    }
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

      int runNo, seqNo;

      if(m_seedFile == "") {
	runNo = m_RunNumber;
	seqNo = m_SequenceNumber;
	++m_SequenceNumber;
      }
      else {

    // index of which event is being generated, only used when reading seeds
    // from a file to generate a set of events
	static int iEvent = 0;

	runNo = m_seeds[iEvent].m_run;
	seqNo = m_seeds[iEvent].m_seqNo;
	++iEvent; 
      }

      // recorde seeds to TDS
      mcevt->initialize(runNo, -1, seqNo);

      // Set run number in EventHeader for consistency with
      // with run number in MCEvent
      SmartDataPtr<Event::EventHeader> header(eventSvc, 
					      EventModel::EventHeader);
      if (header == 0) {
	log << MSG::ERROR << "Error accessing Event Header" << endreq;
      }    

      header->setRun(runNo);

      int multiplier = 1; 
      int dummy = 0; // for 2nd argument to setSeed
      std::map< std::string, HepRandomEngine* >::const_iterator dllEngine;
      for (dllEngine = m_engineMap.begin(); dllEngine != m_engineMap.end(); ++dllEngine ) {
	long theSeed = multiplier * 100000 * ((runNo+1) % 20000) + ((2*seqNo+1) % 100000);
      // note: the following should be DEBUG rather than INFO
      //       but for some reason DEBUG won't print here
      // log << MSG::INFO << "Setting seed for " <<  dllEngine->first << " to " <<  theSeed << endreq;
	dllEngine->second->setSeed(theSeed,dummy);
	++multiplier;
      }
    }
}


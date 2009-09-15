/** @file GlastRandomSvc.cxx 
   @brief Implementation file for GlastRandomSvc

gets adresses 
 of CLHEP random number engines used in Gleam shared libraries
 and sets seeds for them based on run and particle sequence
 number obtained from the MCHeader

 $Header$

 Author: Toby Burnett, Karl Young
*/

#include "GlastRandomSvc.h"
#include "facilities/Util.h"

#include <iterator>
#include <fstream>
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataStoreItem.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/IObjManager.h"
#include "GaudiKernel/IToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"


#include "Event/TopLevel/MCEvent.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include "CLHEP/Random/RandomEngine.h"

// for list of random engines
#include "CLHEP/Random/DualRand.h"
#include "CLHEP/Random/TripleRand.h"
#include "CLHEP/Random/RandEngine.h"
#include "CLHEP/Random/DRand48Engine.h"
#include "CLHEP/Random/Hurd160Engine.h"
#include "CLHEP/Random/Hurd288Engine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Random/RanecuEngine.h"
#include "CLHEP/Random/Ranlux64Engine.h"
#include "CLHEP/Random/RanluxEngine.h"
#include "CLHEP/Random/RanshiEngine.h"
#include "CLHEP/Random/RandGauss.h"


// declare the service factories for the GlastRandomSvc
static SvcFactory<GlastRandomSvc> a_factory;
const ISvcFactory& GlastRandomSvcFactory = a_factory; 

GlastRandomSvc* GlastRandomSvc::instance(){ return s_instance; }
GlastRandomSvc* GlastRandomSvc::s_instance;

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
    declareProperty("RandomEngine",    m_randomEngine="TripleRand");
    declareProperty("RunNumber",       m_RunNumber=-1);
    declareProperty("RunNumberString", m_RunNumberString="$(runName)");
    declareProperty("InitialSequenceNumber", m_InitialSequenceNumber=0);
#if 0 // disable for now: THB
    declareProperty("SeedFile", m_seedFile="");
    declareProperty("EndSeedFile", m_endSeedFile="");
#endif
    declareProperty("autoSeed", m_autoSeed=false); // allow turn off
    s_instance=this; // for access local to GlastSvc
}

GlastRandomSvc::~GlastRandomSvc()  
{
    // this is  repetetive in case finalize was
    // not called.
    finalize();
}


/** @class RanGenFactoryBase
 *
* @brief abstract base class for list of HepRandomEngine factories 
*/
class RanGenFactoryBase {
public:
  virtual ~RanGenFactoryBase(){}
    virtual CLHEP::HepRandomEngine * create() =0;
    virtual std::string name()const =0;
};

/** @class RanGenFactory
 *
* @brief Templated class  to create a particular HepRandomEngine factory 
*/
template <class T> class RanGenFactory : public RanGenFactoryBase{ 
public:
    RanGenFactory(std::string name):m_name(name){}
    virtual CLHEP::HepRandomEngine * create(){return new T;}
    std::string name()const { return m_name;}
private:
    std::string m_name;

};

CLHEP::HepRandomEngine* GlastRandomSvc::createEngine(std::string  engineName) 
{
    // Purpose and Method: Given a string containing a valid CLHEP
    //                     random engine type returns an instance 
    //
    // Inputs:  Name of CLHEP random engine to set
    // Outputs: pointer to a HepRandomEngine, or zero
    // Dependencies: None
    // Restrictions and Caveats:  None

    static RanGenFactoryBase*  factories[] = { 
        new RanGenFactory<CLHEP::DRand48Engine>("DRand48Engine")   ,
        new RanGenFactory<CLHEP::DualRand>("DualRand")             ,
        new RanGenFactory<CLHEP::TripleRand>("TripleRand")         ,
        new RanGenFactory<CLHEP::RandEngine>("RandEngine")         ,
        new RanGenFactory<CLHEP::HepJamesRandom>("HepJamesRandom") ,
        new RanGenFactory<CLHEP::HepJamesRandom>("JamesRandom")    ,
        new RanGenFactory<CLHEP::Hurd160Engine>("Hurd160Engine")   ,
        new RanGenFactory<CLHEP::Hurd288Engine>("Hurd288Engine")   ,
        new RanGenFactory<CLHEP::MTwistEngine>("MTwistEngine")     ,
        new RanGenFactory<CLHEP::RanecuEngine>("RanecuEngine")     ,
        new RanGenFactory<CLHEP::RanshiEngine>("RanshiEngine")     ,
        new RanGenFactory<CLHEP::RanluxEngine>("RanluxEngine")     ,
        new RanGenFactory<CLHEP::Ranlux64Engine>("Ranlux64Engine") 
    };

    for(unsigned int i=0; i< sizeof(factories)/sizeof(void*); ++i){
        if( engineName==factories[i]->name()) return factories[i]->create();
    }
    return 0;
}

StatusCode GlastRandomSvc::initialize () 
{
    // Purpose and Method: initialization: see class description
    //
    // Inputs:  none
    // Outputs: none
    // Dependencies: None
    // Restrictions and Caveats:  None
    StatusCode  status = StatusCode::SUCCESS;

    // Open the message log
    MsgStream log( msgSvc(), name() );

    // Call super-class
    Service::initialize ();

    if( serviceLocator() ) {
        status = serviceLocator()->service("EventDataSvc", m_eventSvc, true );
    }
    if(status.isFailure())
    {
        log << MSG::ERROR << "Could not find EventDataSvc" << endreq;
        return status;
    }


    // Bind all of the properties for this service
    if ( (status = setProperties()).isFailure() ) {
        log << MSG::ERROR << "Failed to set properties" << endreq;
        return status;
    }
    m_SequenceNumber = m_InitialSequenceNumber;

    if( m_RunNumberString!="" ){
        std::string rs(m_RunNumberString.c_str());
        if( (rs.substr(0,2)=="$(") && (rs.substr(rs.size()-1)==")") ) {
            // delimited string should contain an env var -- check it
            const char * rnstring = ::getenv(rs.substr(2, rs.size()-3).c_str());
            if( rnstring !=0) {
                // yes: assume that it is a number
                m_RunNumber = facilities::Util::atoi(rnstring);
                log << MSG::INFO << "Setting run number from environment variable \"" 
                    << rs.substr(2,rs.size()-3) << "\"" <<endreq;
            } 
        }else {
            // string does not define an environment variable
            m_RunNumber = facilities::Util::atoi(rs);
            log << MSG::INFO << "Run number set from numeric string to "
                    << m_RunNumberString << endreq;
        }
    }
    if( m_RunNumber == -1){
        m_RunNumber = 10;
        log << MSG::INFO << "Run number set to default: " << m_RunNumber << endreq;
    } else {
        log << MSG::INFO << "Run number set to: " <<  m_RunNumber << endreq;
    }
    log << MSG::INFO << "==================================================" <<endreq;



    // use the incident service to register begin, end events
    IIncidentSvc* incsvc = 0;
    status = service ("IncidentSvc", incsvc, true);

    if( status.isFailure() ) return status;

    incsvc->addListener(this, "BeginEvent", 0);
    incsvc->addListener(this, "EndEvent", 0);


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

    IToolFactory* toolfactory = 0;

    // Search through all objects (factories?)
    for(IObjManager::ObjIterator it = objManager->objBegin(); it
        !=objManager->objEnd(); ++ it)
    {
        std::string tooltype= (*it)->ident();
        // is it a tool factory?
        const IFactory* factory = objManager->objFactory( tooltype );
        IFactory* fact = const_cast<IFactory*>(factory);
        status = fact->queryInterface( IID_IToolFactory, (void**)&toolfactory );
        if( status.isSuccess() ) 
        {
            std::string fullname = this->name()+"."+tooltype;
            IAlgTool* itool = toolfactory->instantiate(fullname,  this );
            IRandomAccess* ranacc ;
            status =itool->queryInterface(IRandomAccess::interfaceID(), (void**)&ranacc);
            if( status.isSuccess() )
            {
                // Set the Random engine by name
                CLHEP::HepRandomEngine* hr =createEngine(m_randomEngine);
                if( hr==0) 
                {
                    log << MSG::ERROR << "Could not set up the random engine " 
                            << m_randomEngine << endreq;
                    return StatusCode::FAILURE;
                }
                log << MSG::INFO << "Setting CLHEP Engine "<< m_randomEngine
                        << " for " << tooltype << " at " << hr << endreq;
                CLHEP::HepRandomEngine* old = ranacc->setTheEngine(hr);
                // make sure that the old one was not already stored
                for( EngineMap::iterator eit = m_engineMap.begin(); eit != m_engineMap.end(); ++eit)
                {
                    if( eit->second != old) continue;
                    log << MSG::WARNING 
                            << " Previous engine ("<< old << ") for " 
                            << tooltype << " was also set for "<< eit->first << endreq;
                    break;
                }
                // Store its name and address in a map
                m_engineMap[tooltype] = hr;
                IRandomAccess::SetFlag randset = ranacc->getRandSet();
                log << MSG::INFO << "RandGauss setFlag at " << randset << endreq;
                m_setFlagPointers.push_back(randset);
            }

            // Release the tool just allocated
            itool->release();
        }
    }
    // now initialize seeds to be unique for each run 

    log << MSG::DEBUG << "initialize(): calling applySeeds(" << m_RunNumber
        << ", " << m_InitialSequenceNumber << ')' <<endreq;
    if ((m_autoSeed) && (m_RunNumber > 20000)) {
        log << MSG::WARNING << "Should not set run number > 20000 when setting "
            << "seeds event by event.  RunNo is set to " << m_RunNumber
            << " this will result in events that duplicate some run number < "
            << " 20000" << endreq;
    }
    applySeeds(m_RunNumber, m_InitialSequenceNumber);

    return StatusCode::SUCCESS;
}



void GlastRandomSvc::handle(const Incident &inc)
{
     // Purpose and Method: handle "incidents" - i.e. set seeds at beginning of event loop 
    //
    // Inputs:  the incident to handle (begin event)
    // Outputs: new seeds in all the engines
    // Dependencies: None
    // Restrictions and Caveats:  None

    MsgStream log( msgSvc(), name() );

    if( inc.type()=="BeginEvent" ) {

        // See if MCEvent was set up properly
        SmartDataPtr<Event::MCEvent> mcevt(m_eventSvc, EventModel::MC::Event);
        //if (mcevt == 0) {
        if (!mcevt) {
            log << MSG::ERROR << "Error accessing MCEvent" << endreq;
            return;
        }


        int runNo = m_RunNumber, 
            seqNo = m_SequenceNumber++;


        // recorde seeds to TDS
        mcevt->initialize(runNo, -1, seqNo, 0); // last arg is timestamp, will be set in FluxAlg

        // Set run number in EventHeader for consistency with
        // with run number in MCEvent
        SmartDataPtr<Event::EventHeader> header(m_eventSvc, 
            EventModel::EventHeader);
        //if (header == 0) {
        if (!header) {
            log << MSG::ERROR << "Error accessing Event Header: not setting seed" << endreq;
            return;
        }    

        // fill header with run number and event id
        header->setRun(runNo);
        header->setEvent(seqNo);

        if ( m_autoSeed ) {
            log << MSG::DEBUG << "handle(): calling applySeeds(" << runNo
                << ", " << seqNo << ')' << endreq;
            applySeeds(runNo, seqNo);
        }

    }

    if(inc.type()=="EndEvent" && m_endSeedFile.value() != "") {

      static int i = 0;

      m_output << i << std::endl;
      ++i;

      EngineMap::const_iterator dllEngine;
      for (dllEngine = m_engineMap.begin(); dllEngine != m_engineMap.end(); ++dllEngine ) {

    m_output << dllEngine->second->getSeed() << " ";
      }

      m_output << std::endl;
    }

}

void GlastRandomSvc::applySeeds(int runNo, int seqNo)
{
         MsgStream log( msgSvc(), name() );


        int multiplier = 1; 
        int dummy = 0; // for 2nd argument to setSeed
        EngineMap::const_iterator dllEngine;
        for (dllEngine = m_engineMap.begin(); dllEngine != m_engineMap.end(); ++dllEngine ) {
            long theSeed;
            if (m_autoSeed) {
                theSeed = multiplier * 100000 * ((runNo+1) % 20000) + 2*seqNo+1;
            }
            else {
                long bigEnough = 1000000;
                // From Michael Kuss Nov 30, 2007
                theSeed = 2 * ( (runNo+1) % bigEnough + multiplier * bigEnough * (seqNo+1) ) + 1;
            }
                

            log << MSG::DEBUG;
            if(log.isActive() ){
                log << "Setting seed for " <<  dllEngine->first << " to " <<  theSeed ;
            }
            log << endreq;
            
            dllEngine->second->setSeed(theSeed,dummy);
            ++multiplier;
        }
        // Gaussian are always produced in pairs, one of them being returned,
        // the other one being cached for a second call.  Setting the seed
        // doesn't flush the cache!
        CLHEP::RandGauss::setFlag(false);

        // now make sure the external ones (if any) are done too!
        for( std::vector< IRandomAccess::SetFlag>::iterator it = m_setFlagPointers.begin();
            it != m_setFlagPointers.end(); ++it)
        {
            (**it)(false);
        }

}

StatusCode GlastRandomSvc::finalize () 
{
    // Purpose and Method: finish up, delete the list of engines
    //
    // Inputs:  none
    // Outputs: none
    // Dependencies: None
    // Restrictions and Caveats:  None

    for (EngineMap::iterator dllEngine = m_engineMap.begin(); 
        dllEngine != m_engineMap.end(); ++dllEngine ) {
        delete dllEngine->second;       
    }
    m_engineMap.clear();
    return StatusCode::SUCCESS;
}


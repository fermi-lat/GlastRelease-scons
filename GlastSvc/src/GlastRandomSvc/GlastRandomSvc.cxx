/*
 *
 * @file GlastRandomSvc.cxx 
 *
 * @brief Implementation file for GlastRandomSvc
 *
 * gets adresses 
 * of CLHEP random number engines used in Gleam shared libraries
 * and sets seeds for them based on run and particle sequence
 * number obtained from the MCHeader
 * 
 * $Header$
 *
 * Author: Toby Burnett, Karl Young, Michael Kuss
 *
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

GlastRandomSvc::GlastRandomSvc(const std::string& name,ISvcLocator* svc)
    : Service(name,svc)
{
    // Purpose and Method: Constructor - Declares and sets default properties
    //                     
    // Inputs: service name and locator 
    //         
    // Outputs: None
    // Dependencies: None
    // Restrictions and Caveats:  None

    // declare the properties
    declareProperty("RandomEngine",          m_randomEngine="TripleRand");
    declareProperty("RunNumber",             m_runNumber=10);
    declareProperty("RunNumberString",       m_runNumberString="$(runName)");
    declareProperty("InitialSequenceNumber", m_initialSequenceNumber=0);
    declareProperty("NumberOfRuns",          m_numberOfRuns=0);
    declareProperty("Seed",                  m_seed=0);
    //    declareProperty("eventSeeded",        m_eventSeeded=false);
    //    declareProperty("autoSeed",              m_autoSeed=false);
    declareProperty("autoSeed",              m_eventSeeded=false);
    declareProperty("pedantic",              m_pedantic=true);
#if 0 // disable for now: THB
    declareProperty("SeedFile", m_seedFile="");
    declareProperty("EndSeedFile", m_endSeedFile="");
#endif
    s_instance=this; // for access local to GlastSvc
}

GlastRandomSvc::~GlastRandomSvc()  
{
    // this is  repetetive in case finalize was
    // not called.
    finalize();
}


/*
 * @class RanGenFactoryBase
 *
 * @brief abstract base class for list of HepRandomEngine factories 
 */
class RanGenFactoryBase {
public:
  virtual ~RanGenFactoryBase(){}
    virtual CLHEP::HepRandomEngine * create() =0;
    virtual std::string name()const =0;
};

/*
 * @class RanGenFactory
 *
 * @brief Templated class to create a particular HepRandomEngine factory 
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

    // Call super-class, here the job options parameters obtain their values
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
    m_sequenceNumber = m_initialSequenceNumber;

    // initialize m_numberOfRuns here, conditional initialization of Gaudi
    // properties is not possible (in the constructor)
    if ( m_numberOfRuns == 0 )
        m_numberOfRuns = m_eventSeeded ? 20000 : 20000000;
    
    if( m_runNumberString!="" ){
        std::string rs(m_runNumberString.c_str());
        if( (rs.substr(0,2)=="$(") && (rs.substr(rs.size()-1)==")") ) {
            // delimited string should contain an env var -- check it
            const char * rnstring = ::getenv(rs.substr(2, rs.size()-3).c_str());
            if( rnstring !=0) {
                // yes: assume that it is a number
                m_runNumber = facilities::Util::atoi(rnstring);
                log << MSG::INFO
                    << "Setting run number from environment variable \"" 
                    << rs.substr(2,rs.size()-3) << "\"" <<endreq;
            } 
        }else {
            // string does not define an environment variable
            m_runNumber = facilities::Util::atoi(rs);
            log << MSG::INFO << "Run number set from numeric string to "
                    << m_runNumberString << endreq;
        }
    }

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
    status = serviceLocator()->service("ApplicationMgr", objManager);
    if ( status.isFailure() ) {
        log << MSG::ERROR << "Unable to locate ObjectManager Service" << endreq;
        return status;
    }

    IToolFactory* toolfactory = 0;

    m_engineNum = 0;

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
            status =itool->queryInterface(IRandomAccess::interfaceID(),
                                          (void**)&ranacc);
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
                ++m_engineNum;
                CLHEP::HepRandomEngine* old = ranacc->setTheEngine(hr);
                // make sure that the old one was not already stored
                for( EngineMap::iterator eit = m_engineMap.begin();
                     eit != m_engineMap.end(); ++eit)
                {
                    if( eit->second != old) continue;
                    // MWK: I downgrade this from WARNING to INFO.
                    // We know that linux/gcc instantiates only one engine,
                    // and the fact is handled in the code.
                    log << MSG::INFO 
                        << " Previous engine ("<< old << ") for " << tooltype
                        << " was also set for "<< eit->first << endreq;
                    --m_engineNum;
                    break;
                }
                // Store its name and address in a map
                m_engineMap[tooltype] = hr;
                IRandomAccess::SetFlag randset = ranacc->getRandSet();
                log << MSG::INFO << "RandGauss setFlag at " << randset <<endreq;
                m_setFlagPointers.push_back(randset);
            }

            // Release the tool just allocated
            itool->release();
        }
    }

    log << MSG::INFO << "number of engines: " << m_engineMap.size()
        << " (unique instantiations: " << m_engineNum << ")" << endreq;

    log << MSG::INFO << "autoSeed set to: " << (m_eventSeeded?"true":"false");
    if ( log.isActive() )
        if ( m_eventSeeded )
            log << " (unique seed for each event)";
        else
            log << " (one seed at the beginning of the run)";
    log << endreq;

    log << MSG::INFO << "RandomEngine set to: " << m_randomEngine << endreq;
    log << MSG::INFO << "RunNumber set to: " << m_runNumber << endreq;
    log << MSG::INFO << "InitialSequenceNumber set to: "
        << m_initialSequenceNumber << endreq;
    log << MSG::INFO << "NumberOfRuns set to: " << m_numberOfRuns<<endreq;

    // Setting the parameter "Seed" is invalid:
    // 1) for event seeded runs
    // 2) for runs with multiple engine instantiations
    if ( m_seed > 0 ) {
        if ( m_eventSeeded || m_engineNum != 1 ) {
            log << MSG::ERROR << "Explicitely setting the seed (here to "
                << m_seed << ") in the job options doesn't make sense:";
            //            m_seed = 0;
            if ( log.isActive() ) {
                if ( m_eventSeeded )
                    log << " *) for runs with unique seed for each event";
                if ( m_engineNum != 1 )
                    log << " *) for more than 1 engine instantiation";
            }
            log << endreq;
            return StatusCode::FAILURE;
        }
        else
            log << MSG::INFO << "Seed set to: " << m_seed << endreq;
    }

    // pick up ApplicationMgr.EvtMax
    IProperty* propMgr = 0;
    status = serviceLocator()->service("ApplicationMgr", propMgr);
    if ( status.isFailure() ) {
        log << MSG::ERROR << "Unable to locate PropertyManager Service"<<endreq;
        return status;
    }
    IntegerProperty AppMgrEvtMax("EvtMax", 0);
    status = propMgr->getProperty(&AppMgrEvtMax);
    if ( status.isFailure() ) {
        log << MSG::ERROR << "could not retrieve ApplicationMgr.EvtMax"<<endreq;
        return status;
    }
    log << MSG::INFO << "ApplicationMgr.EvtMax set to: " <<AppMgrEvtMax<<endreq;

    const int twoTo30 = 1073741824;
    //        const int twoTo30 = 1000000000;
    // subtract 1 from m_maxEventId: e.g. 8 numbers in the range 0
    // to 7
    m_maxEventId = twoTo30 / ( m_numberOfRuns * m_engineNum ) - 1;
    log << MSG::INFO;
    if ( log.isActive() ) {
        log << "maximum allowed event id";
        if ( ! m_eventSeeded )
            log << " (InitialSequenceNumber)";
        log << ": " << m_maxEventId;
    }
    log << endreq;

    // sanity checks
    bool exitFlag = false;

    // RunNumber should not be negative
    if ( m_runNumber < 0 ) {
        log << MSG::WARNING << "chosen RunNumber is " << m_runNumber
            << ": a negative RunNumber is insane" << endreq;
        exitFlag = true;
    }

    // InitialSequenceNumber should not be negative
    if ( m_initialSequenceNumber < 0 ) {
        log << MSG::WARNING << "chosen InitialSequenceNumber is "
            << m_initialSequenceNumber
            << ": a negative InitialSequenceNumber is insane" << endreq;
        exitFlag = true;
    }

    // RunNumber should be smaller than NumberOfRuns
    if ( m_runNumber >= m_numberOfRuns ) {
        log << MSG::WARNING << "chosen runNumber " << m_runNumber
            << " is not smaller than the total number of runs "
            << m_numberOfRuns << endreq;
        exitFlag = true;
    }

    // the maximum event id in a run (given by InitialSequenceNumber+
    // ApplicationMgr.EvtMax-1 for event seeded runs, given by
    // InitialSequenceNumber for run seeded runs) should not exceed
    // m_maxEventId.
    const int maxEventIDofRun = m_initialSequenceNumber
        + ( m_eventSeeded ? AppMgrEvtMax - 1 : 0 );
    log << MSG::DEBUG << "maximum event id that will be used for seeding in "
        << "this run: " << maxEventIDofRun << endreq;
    if ( maxEventIDofRun > m_maxEventId ) {
        log << MSG::WARNING << "Chosen NumberOfRuns and the determined number "
            << "of engine instantiations allow for a maximum event id of "
            << m_maxEventId;
        if ( log.isActive() ) {
            if ( m_eventSeeded )
                log <<", but InitialSequenceNumber+ApplicationMgr.EvtMax-1 is ";
            else
                log << ", but InitialSequenceNumber is ";
            log << maxEventIDofRun;
        }
        log << endreq;
        exitFlag = true;
    }            

    // check for consistency in engine numbers.
    // The number of engine instantiations should be either the number of
    // engines (i.e. Windows/VisualC++) or 1 (Linux/gcc).
    if ( m_engineNum != 1 && m_engineNum != m_engineMap.size() ) {
        log << MSG::WARNING << "The number of engine instantiations, here "
            << m_engineNum << ", should be equal to either 1 or the total "
            << "number of engines, here " << m_engineMap.size()
            << ".  Everything else is weird!" << endreq;
        exitFlag = true;
    }

    if ( exitFlag ) {
        log << MSG::WARNING << "The preceeding warnings highlight settings "
            << "which may cause duplication of events." << endreq;
        if ( m_pedantic ) {
            log << MSG::ERROR << "If you don't care for duplicated events, or "
                << "have taken other measures to avoid it, you can inhibit "
                << "termination by specifying 'GlastRandomSvc.pedantic=false;' "
                << "in your jobOptions file." << endreq;
            log << MSG::ERROR << "GlastRandomSvc will exit with "
                << "StatusCode::FAILURE now!" << endreq;
            return StatusCode::FAILURE;
        }
        else
            log << MSG::WARNING<<"Termination disabled by setting pedantic to: "
                << (m_pedantic?"true":"false")
                << ".  You are on your own now!" << endreq;
    }

    log << MSG::DEBUG << "initialize(): calling applySeeds(" << m_runNumber
        << ", " << m_initialSequenceNumber << ')' <<endreq;
    applySeeds(m_runNumber, m_initialSequenceNumber);

    return StatusCode::SUCCESS;
}



void GlastRandomSvc::handle(const Incident &inc)
{
    // Purpose and Method: handle "incidents" - i.e. set seeds at beginning of
    // event loop 
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

        const int runNo = m_runNumber;
        const int seqNo = m_sequenceNumber++;

        // record seeds to TDS
        mcevt->initialize(runNo, -1, seqNo, 0);
        // last arg is timestamp, will be set in FluxAlg

        // Set run number in EventHeader for consistency with
        // run number in MCEvent
        SmartDataPtr<Event::EventHeader> header(m_eventSvc, 
            EventModel::EventHeader);
        //if (header == 0) {
        if (!header) {
            log << MSG::ERROR
                << "Error accessing Event Header: not setting seed" << endreq;
            return;
        }    

        // fill header with run number and event id
        header->setRun(runNo);
        header->setEvent(seqNo);

        if ( m_eventSeeded ) {
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
        for ( dllEngine=m_engineMap.begin(); dllEngine!=m_engineMap.end();
              ++dllEngine )
            m_output << dllEngine->second->getSeed() << " ";
        m_output << std::endl;
    }
}

void GlastRandomSvc::applySeeds(int runNo, int seqNo) {
    MsgStream log(msgSvc(), name());

    int engineNum = 1;
    const int luxury = 0; // 2nd argument to setSeed
    EngineMap::const_iterator dllEngine;
    for ( dllEngine=m_engineMap.begin(); dllEngine!=m_engineMap.end();
          ++dllEngine ) {
        long theSeed;
        if ( m_seed )
            theSeed = m_seed;
        else
            theSeed = 2 * (
                           (runNo%m_numberOfRuns)*m_engineNum*(m_maxEventId+1)
                           + (engineNum-1)*(m_maxEventId+1) + seqNo
                           )
                + 1;
        log << MSG::DEBUG << "Setting seed for " <<  dllEngine->first << " to "
            <<  theSeed << endreq;
        dllEngine->second->setSeed(theSeed, luxury);
        // use engineNum only if there are multiple engine instantiations
        // (Windows/VisualC++)
        if ( m_engineNum > 1 )
            ++engineNum;
    }
    log << MSG::DEBUG;
    if ( log.isActive() ) {
        log << "reading back the seeds:";
        for ( dllEngine=m_engineMap.begin(); dllEngine!=m_engineMap.end();
              ++dllEngine )
            log << MSG::DEBUG << std::endl << "    seed of "
                <<  dllEngine->first << ": " << dllEngine->second->getSeed();
    }
    log << endreq;

    // Gaussian are always produced in pairs, one of them being returned,
    // the other one being cached for a second call.  Setting the seed
    // doesn't flush the cache!
    CLHEP::RandGauss::setFlag(false);

    // now make sure the external ones (if any) are done too!
    for ( std::vector<IRandomAccess::SetFlag>::iterator
              it=m_setFlagPointers.begin();
          it!=m_setFlagPointers.end(); ++it ) {
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

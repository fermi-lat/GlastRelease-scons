/** 
* @file FluxSvc.cxx
* @brief definition of the class FluxSvc
*
*  $Header$
*  Original author: Toby Burnett tburnett@u.washington.edu
*/

#include "FluxSvc/IRegisterSource.h"

#include "facilities/Timestamp.h"
#include "facilities/Observer.h"

#include "celestialSources/SpectrumFactoryLoader.h"


#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/IObjManager.h"
#include "GaudiKernel/IToolFactory.h"
#include "GaudiKernel/IAlgManager.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/IParticlePropertySvc.h"

#include "CLHEP/Random/Random.h"

#include "flux/Flux.h"
#include "flux/FluxMgr.h"
#include "flux/rootplot.h"
#include "flux/ISpectrumFactory.h"

#include <algorithm>
#include <iterator>
#include <sstream>
#include <iomanip>
using astro::GPS;
/** 
* \class FluxSvc
*
* \brief Service that implements the IFluxSvc interface, to return an IFlux object.
*  FluxSvc handles the creation and interfacing with Flux objects.  
* \author Toby Burnett tburnett@u.washington.edu
* 
* $Header$
*/

// includes
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IRunable.h"
#include "GaudiKernel/Property.h"

#include "FluxSvc/IFluxSvc.h"
#include <list>

//forward declarations
template <class TYPE> class SvcFactory;
class IFlux;  // interface
class FluxMgr;  // actual manager
class IParticlePropertySvc; 
class IAppMgrUI;

class FluxSvc : 
    virtual public Service, 
    virtual public IFluxSvc,
    virtual public IRunable
{  
public:

    /// return pointer to a flux object
    StatusCode source(std::string name, IFlux*&);

    /// set pointer to a flux object, constructed from set of names
    StatusCode compositeSource(std::vector<std::string> names, IFlux*& flux);

    /// return a list of possible names
    std::list<std::string> fluxNames()const;

    /// add a new SpectrumFactory
    virtual void addFactory(std::string name, const ISpectrumFactory* factory );

    /// pass a specific amount of time
    virtual void pass ( double t);


    /// return pointer to the random engine that FluxSvc uses
    virtual HepRandomEngine* getRandomEngine();
#if 0
    /// create a set of display windows using rootplot.
    void rootDisplay(std::vector<const char*> arguments);
#endif
    virtual void rootDisplay(std::vector<std::string> arguments);;

    /// attach an external observer to GPS
    void attachGpsObserver(Observer* anObserver);

    ///return the pointer to the current IFlux object
    IFlux* currentFlux();

    /// name of the flux
    std::string fluxName()const;

    /// set the glast tilt angles for explicit, static rocking
    /// the angles correspond to a rotation about the x axis followed by the z.
    void setExplicitRockingAngles(double ang1, double ang2);

    /// get the angular values of the satellite
    std::pair<double,double> getExplicitRockingAngles();

    ///this transforms glast-local (cartesian) vectors into galactic (cartesian) vectors
    HepRotation transformGlastToGalactic(double time)const;

    HepRotation transformToGlast(double seconds,GPS::CoordSystem index)const;

    /// get the current satellite location
    std::pair<double,double> location();

    /// return pointer to our GPS instance
    GPS* GPSinstance();

    /// return a string which uniquely identifies the source
    std::string uniqueIDString()const;

    /// Set the satellite rocking mode:
    ///0=NONE, 1=UPDOWN(up in the northern hemisphere, down in the southern,
    ///2=SLEWING(like updown, but not discontinuous at the equator),
    ///3=ONEPERORBIT (rock norh one orbit, south the next,
    ///4=EXPLICIT (use the internal rotangles rotation angles (this should be set through setOrientation)).
    ///5 = POINT:  Explicit pointing direction given - setExplicitRockingAngles are (l,b).
    ///6 = HISTORY - Filename given to stand for a pre-recorded pointing history.  Use the setPointingHistoryFile function.
    std::vector<double> setRockType(int rockType, double rockAngle = 35.);

    /// set the desired pointing history file to use:
    void setPointingHistoryFile(std::string fileName);

    ///this should return the source file names, along with the contained sources.
    std::vector<std::pair< std::string ,std::list<std::string> > > sourceOriginList() const;


    bool insideSAA() { return m_insideSAA;}
    /// for the IRunnable interfce
    virtual StatusCode run();

    //------------------------------------------------------------------
    //  stuff required by a Service

    /// perform initializations for this service. 
    virtual StatusCode initialize ();

    /// perform the finalization, as required for a service.
    virtual StatusCode finalize ();

    /// Query interface
    virtual StatusCode queryInterface( const IID& riid, void** ppvUnknown );

protected: 

    /// Standard Constructor
    FluxSvc ( const std::string& name, ISvcLocator* al );

    /// destructor
    virtual ~FluxSvc ();

private:
    //! 

    IParticlePropertySvc* m_partSvc;

    /// Allow SvcFactory to instantiate the service.
    friend class SvcFactory<FluxSvc>;

    FluxMgr * m_fluxMgr;
    /// the user-defined list of acceptable XML sources (from JobOptions.txt)
    //std::string m_source_lib;
    std::vector<std::string> m_source_lib;
    /// the default XML file name (from JobOptions.txt)
    std::string m_source_lib_default;
    /// set dtd to use.
    std::string m_dtd_file;
    /// the "current" flux object
    IFlux* m_currentFlux;

    /// Reference to application manager UI
    IAppMgrUI*    m_appMgrUI;
    IntegerProperty m_evtMax;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /** @class Times
        @brief encapsulate the various times of interest
        */
    class Times {
    public:
        void initialize(MsgStream& log){
            // the launch: if set, all times will be added to it
            
            m_start = convert(m_startDate);
            m_launch= convert(m_launchDate);
            if( m_launch==0 ) m_launch=m_start;

            // now add StartTime as offset to either the StartDate or LaunchDate
            if(! m_startTimeEnvVar.value().empty()) {
                const char* envar=::getenv(m_startTimeEnvVar.value().c_str());
                if( envar==0 ){
                    log << MSG::WARNING << "Env var " << m_startTimeEnvVar.value() 
                        << " requested for start time, not found" << endreq;
                } else {
                    m_start = ::atof( envar );
                    log << MSG::INFO << "Setting start time from environment variable " 
                        << m_startTimeEnvVar.value() << " to "
                        << m_start << endreq; 
                }
            } else {
                m_start += m_startTime.value();
            }

            if( m_deltaTime>0 && m_endTime==0 )  m_endTime=m_start+m_deltaTime;
            // set the basic time here: it will be incremented by the flux object
            GPS::instance()->time(m_start);
            m_end = m_endTime;
            if( m_deltaTime>0) m_end=m_start+m_deltaTime;

        }
        double convert(std::string date){
            using astro::JulianDate;
            if(date.empty())return 0;
            // parse a string of the form "2004-09-03 18:00" using the Timestamp class.
            facilities::Timestamp jt(date);
            return (JulianDate(jt.getJulian())-JulianDate::missionStart())*JulianDate::secondsPerDay;
        }
        // properties get set by gaudi 
        StringProperty m_launchDate; ///
        StringProperty m_startDate;
        DoubleProperty m_startTime;
        DoubleProperty m_endTime;
        DoubleProperty m_deltaTime;
        StringProperty m_startTimeEnvVar;

        double m_launch;
        double m_start;
        double m_end;
        double m_current;
        double launch(){return m_launch;}
        double start(){return m_start;}
        double current(){ return GPS::instance()->time();}
        double end(){return m_end;}
        double offset(double t){return t-m_launch;}
    } m_times;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ObserverAdapter< FluxSvc > m_observer; //observer tag
    int askGPS(); // the function that will be called
    bool m_insideSAA; 

    DoubleProperty m_expansionFactor;


};

// declare the service factories for the FluxSvc
static SvcFactory<FluxSvc> a_factory;
const ISvcFactory& FluxSvcFactory = a_factory;


static std::string default_source_library("$(FLUXROOT)/xml/source_library.xml");
static std::string default_dtd_file("$(FLUXROOT)/xml/source.dtd");

// ------------------------------------------------
// Implementation of the FluxSvc class
// ------------------------------------------------
/// Standard Constructor
FluxSvc::FluxSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc), m_currentFlux(0), m_insideSAA(false)
, m_expansionFactor(1.0)
{

    declareProperty("source_lib" , m_source_lib); 
    declareProperty("dtd_file"   , m_dtd_file=default_dtd_file);
    declareProperty("EvtMax"     , m_evtMax=0);

    declareProperty("StartTime"   , m_times.m_startTime=0);
    declareProperty("EndTime"     , m_times.m_endTime=0);
    declareProperty("DeltaTime"   , m_times.m_deltaTime=0);
    declareProperty("StartDate"   , m_times.m_startDate="");
    declareProperty("LaunchDate"  , m_times.m_launchDate="");
    declareProperty("StartTimeEnvVar", m_times.m_startTimeEnvVar="");


#if 0 // disable this for now, it is not consistent with CompositeSource
    declareProperty("ExpansionFactor"   , m_expansionFactor=1.0);
#endif

}

std::list<std::string> FluxSvc::fluxNames()const{
    return m_fluxMgr->sourceList();
}

StatusCode FluxSvc::source(std::string name, IFlux*& flux) {
    std::list<std::string> source_list( fluxNames() );
    std::list<std::string> source_list2( SpectrumFactoryTable::instance()->spectrumList() );

    if( std::find(source_list.begin(), source_list.end(), name) == source_list.end() 
        &&(std::find(source_list2.begin(), source_list2.end(), name) == source_list2.end()))
        return StatusCode::FAILURE;
    flux =  new Flux(name);
    m_currentFlux = flux;    
    return StatusCode::SUCCESS;
}

StatusCode FluxSvc::compositeSource(std::vector<std::string> names, IFlux*& flux) {
    flux =  new Flux(names);
    if( flux->currentEvent()==0) return StatusCode::FAILURE;
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

    m_times.initialize(log);

    status = serviceLocator()->queryInterface(IID_IAppMgrUI, (void**)&m_appMgrUI);

    // If source library was not set, put in default
    if( m_source_lib.empty() ){
        m_source_lib.push_back(default_source_library);
        log << MSG::INFO << "Set source library list to " << default_source_library << endreq;
    }

    try {
        // create a FluxMgr object which will then be available.
        m_fluxMgr = new FluxMgr(m_source_lib, m_dtd_file);
    }catch(...){
        return StatusCode::FAILURE;
    }

    Flux::mgr(m_fluxMgr); // tell our Flux object

    // check that it was made properly
    if( m_fluxMgr->sourceList().empty()) {
        log << MSG::ERROR  << "Did not initialize properly: no sources detected" << endreq;
        status = StatusCode::FAILURE;
    }

    m_fluxMgr->setExpansion(m_expansionFactor);

    if ( service("ParticlePropertySvc", m_partSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the ParticlePropertySvc!" << endreq;
        return StatusCode::FAILURE;
    }

    log << MSG::INFO << "Registering factories external to flux: ";
    SpectrumFactoryLoader externals;
    std::vector<std::string> flux_names(externals.names());

    std::copy( flux_names.begin(), flux_names.end(), 
        std::ostream_iterator<std::string>(log.stream(), ", "));
    log  << endreq;



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

    IToolFactory* toolfactory = 0;
    
    // search throught all objects (factories?)
    for(IObjManager::ObjIterator it = objManager->objBegin(); it !=objManager->objEnd(); ++ it){

        std::string tooltype= (*it)->ident();
        // is it a tool factory?
        const IFactory* factory = objManager->objFactory( tooltype );
        IFactory* fact = const_cast<IFactory*>(factory);
        status = fact->queryInterface( IID_IToolFactory, (void**)&toolfactory );
        if( status.isSuccess() ) {

            IAlgTool* itool = toolfactory->instantiate(name()+"."+tooltype,  this );
            IRegisterSource* ireg;
            status =itool->queryInterface( IRegisterSource::interfaceID(), (void**)&ireg);
            if( status.isSuccess() ){
                log << MSG::INFO << "Registering sources in " << tooltype << endreq;
                ireg->registerMe(this);
            }
            log << MSG::DEBUG << "Releasing the tool " << tooltype << endreq;
            itool->release();
        }

    }

    // attach an observer to be notified when orbital position changes
    // set callback to be notified when the position changes
    m_observer.setAdapter( new ActionAdapter<FluxSvc>
        (this, &FluxSvc::askGPS) );

    attachGpsObserver(&m_observer);


    return StatusCode::SUCCESS;
}
//------------------------------------------------------------------------
int FluxSvc::askGPS()
{
    astro::EarthCoordinate pos = GPS::instance()->earthpos();
    bool inside = pos.insideSAA();

    if( m_insideSAA == inside) return 0; // no change
    MsgStream log( msgSvc(), name() );
    log << MSG::INFO
        << (!inside? " leaving" : "entering") 
        << " SAA at "  << GPS::instance()->time() << endreq;
    m_insideSAA = inside;

    return 0; // can't be void in observer pattern
}


/// return pointer to the random engine
HepRandomEngine* FluxSvc::getRandomEngine(){
    return  HepRandom::getTheEngine();
};

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
#if 0
void FluxSvc::rootDisplay(std::vector<const char*> arguments){
    rootplot abc(arguments, m_fluxMgr);
}
#endif

void FluxSvc::rootDisplay(std::vector<std::string> arguments){
    rootplot abc(arguments, m_fluxMgr);   
}


void FluxSvc::attachGpsObserver(Observer* anObserver)
{
    GPS::instance()->notification().attach( anObserver );
}


/// return pointer to the GPS instance of FluxSVc
GPS* FluxSvc::GPSinstance(){ return GPS::instance();}


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
    std::stringstream t;
    t << m_currentFlux->numSource();
    return m_currentFlux->name() + t.str();
}


void FluxSvc::setExplicitRockingAngles(double ang1, double ang2){
    m_fluxMgr->setExplicitRockingAngles(std::make_pair<double,double>(ang1,ang2));
}

/// get the angular values of the satellite
std::pair<double,double> FluxSvc::getExplicitRockingAngles(){
    return m_fluxMgr->getExplicitRockingAngles();
}

void FluxSvc::setPointingHistoryFile(std::string fileName){
    m_fluxMgr->setPointingHistoryFile(fileName);
}

HepRotation FluxSvc::transformToGlast(double seconds,GPS::CoordSystem index)const{
    return m_fluxMgr->transformToGlast(seconds,index);
}

HepRotation FluxSvc::transformGlastToGalactic(double time)const{
    return m_fluxMgr->transformGlastToGalactic(time);
}

std::pair<double,double> FluxSvc::location(){
    return m_fluxMgr->location();
}
/// this sets the rocking mode in GPS.
std::vector<double> FluxSvc::setRockType(int rockType, double rockAngle){
    return m_fluxMgr->setRockType(rockType, rockAngle);
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
    int eventNumber= 0;
    
    // access the starting time from the job properties
    { bool noend=true;
    log << MSG::INFO << "Runable interface starting event loop as :" ; 
    if( m_evtMax>0)  { log << " MaxEvt = " << m_evtMax; noend=false;  }
    if( m_times.start()>0) { log << " StartTime= launch+" << m_times.start() -m_times.launch(); }
    if( m_times.end()>0  ) { log << " EndTime=launch+ " << m_times.end()-m_times.launch(); noend=false; }
    log << endreq;

    if(noend) { 
        log << MSG::ERROR<< "No end condition specified: will not process any events!" << endreq; 
        return StatusCode::FAILURE;
    }
    }
    if( m_times.end()>0 && m_times.start() > m_times.end()){
        log << MSG::ERROR << "Start time after end time!" << endreq;
        return StatusCode::FAILURE;
    }
    int last_fraction=0;
    // loop: will quit if either limit is set, and exceeded
    while( (m_evtMax==0  || m_evtMax>0 &&  eventNumber < m_evtMax)
        && (m_times.end()==0 || m_times.end()>0 && m_times.current() < m_times.end()) ) {

            double efrac =   (m_evtMax>0? 100.*eventNumber/m_evtMax: 0.0), 
                tfrac =   (m_times.end()>0? 100.*(m_times.current()-m_times.start())/(m_times.end()-m_times.start()) : 0.0) ;

            int percent_complete= static_cast<int>(  std::max( efrac, tfrac)  );
            if( percent_complete!=last_fraction){
                last_fraction=percent_complete;
                if( percent_complete<10 || percent_complete%10 ==0){
                    log << MSG::INFO <<   percent_complete << "% complete: "
                        << " event "<< eventNumber<<",  time=launch+ "<< (m_times.current()-m_times.launch()) << endreq;
                }

            }

            status =  m_appMgrUI->nextEvent(1); // currently, always success

            // the single event may have created a failure. Check the ErrorCount propery of the Top alg.
            if( theAlgorithm !=0) theAlgorithm->getProperty(&errorProperty);
            if( status.isFailure() || errorProperty.value() > 0){
                status = StatusCode::FAILURE;
            }

            if( status.isFailure()) break;
            eventNumber ++;
        }
        if( status.isFailure()){
            log << MSG::ERROR << "Terminating FluxSvc loop due to error" << endreq;

        }else if( m_times.end()>0 && m_times.current() >= m_times.end() ) {
            log << MSG::INFO << "Loop terminated by time " << endreq;
        }else {
            log << MSG::INFO << "Processing loop terminated by event count" << endreq;
        }
        log << MSG::INFO << "End after "<< eventNumber << " events, time = " << m_times.current() << endreq;
        return status;
    }


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
#include "facilities/Util.h"

#include "celestialSources/SpectrumFactoryLoader.h"


#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/IAlgManager.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/IParticlePropertySvc.h"

#include "FluxObs.h"

#include "CLHEP/Random/Random.h"

#include "GlastSvc/GlastRandomSvc/IRandomAccess.h"

#include "astro/SkyDir.h"
#include "astro/EarthOrbit.h"
#include "astro/EarthCoordinate.h"

#include "flux/Flux.h"
#include "flux/FluxMgr.h"
#include "flux/rootplot.h"
#include "flux/ISpectrumFactory.h"
#include "flux/Spectrum.h"

#include <cassert>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <iostream>

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
    virtual CLHEP::HepRandomEngine* getRandomEngine();
    virtual void rootDisplay(std::vector<std::string> arguments);;

    /// attach an external observer to GPS
    void attachGpsObserver(Observer* anObserver);

    ///return the pointer to the current IFlux object
    IFlux* currentFlux();

    /// name of the flux
    std::string fluxName()const;

    /// set the pointing direction 
    void setPointingDirection(const astro::SkyDir& dir);

    /// get the angular values of the satellite
    std::pair<double,double> getExplicitRockingAngles();

    ///this transforms glast-local (cartesian) vectors into galactic (cartesian) vectors
    CLHEP::HepRotation transformGlastToGalactic(double time)const;

    CLHEP::HepRotation transformToGlast(double seconds,GPS::CoordSystem index)const;

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
    std::vector<double> setRockType(astro::GPS::RockType rockType, double rockAngle = 35.);

    ///this should return the source file names, along with the contained sources.
    std::vector<std::pair< std::string ,std::list<std::string> > > sourceOriginList() const;


    bool insideSAA() { return m_insideSAA;}


    /// for the IRunnable interfce
    virtual StatusCode run();

    double endruntime(); ///< access the end of run time 

    virtual void setFilterCone(std::vector<double> cone); ///< set filter cone parameters (ra, dec, radius)

    /// set aligmnment for Glast. 
    /// @param qx, qy, qz x, y, z rotation angles (radians) asssumed small
    /// @param misalign [false] set true to apply as a misalignment
    virtual void setAlignmentRotation(double qx, double qy, double qz, bool misalign);

    /// set the SAA boundary
    /// @param boundard set of (lat, lon) pairs of the polygon
    virtual void setSAABoundary(const std::vector<std::pair<double, double> > & boundary);


    //------------------------------------------------------------------
    //  stuff required by a Service

    /// perform initializations for this service. 
    virtual StatusCode initialize ();

    /// perform the finalization, as required for a service.
    virtual StatusCode finalize ();

    /// Query interface
    virtual StatusCode queryInterface( const InterfaceID& riid, void** ppvUnknown );

protected: 

    /// Standard Constructor
    FluxSvc ( const std::string& name, ISvcLocator* al );

    /// destructor
    virtual ~FluxSvc ();

private:
    //! 

    IParticlePropertySvc* m_partSvc;

    IRandomAccess *m_randTool;

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

    FluxObs *m_fluxObs;
    IToolSvc *m_toolSvc; // to handle observer


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
            double offset = m_startTime.value();
            double delta  = m_deltaTime.value();
	    int comma = -1;
            if(! m_startTimeEnvVar.value().empty()) {
	      //user provided info for this parameter, so...
                std::string ev(m_startTimeEnvVar.value());
		try {
		  //... first try the old way, needed for gaudi versions 
		  // where it is not expanded on the fly when calling .value()
		  int success = facilities::Util::expandEnvVar(&ev);
		  comma =  ev.find(',');
		} catch(facilities::Untranslatable ex)
		  {
		    //assume that new gaudi expanded the env var already
		    //In that case we require here that a comma be present
		    comma =  ev.find(',');
		    if(comma==-1){
		      log<<MSG::FATAL<<"Incorrect value \""<<ev<<"\" for property FluxSvc.startTimeEnvVar : it should be an environment variable or a pair of numbers separated by a comma"<<endreq; 
		      throw std::runtime_error("");
		    }
		  }
		offset = ::atof( ev.substr(0,comma>0?comma+1:-1).c_str() );
		log << MSG::INFO << "Setting start time offset from environment variable " 
		    << m_startTimeEnvVar << " to "  << offset; 
		if( comma>0){ 
		  delta = ::atof(ev.substr(comma+1).c_str());
		  log << "\n\t\t\tsetting deltat to " << delta; 
		}
		log << endreq;
	    }
            m_start += offset;
            if( delta >0 && m_endTime==0 )  m_endTime=m_start+delta;
            // set the basic time here: it will be incremented by the flux object
            GPS::instance()->time(m_start);
            m_end = m_endTime;
            if( m_deltaTime>0) m_end=m_start+delta;

	    log << MSG::INFO << "init: start time = " << std::setprecision(12)
		<< m_start << " sec, end time = " << m_end << " sec, delta = " 
		<< m_end-m_start << " sec" << endreq;

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
    DoubleProperty m_sampleInterval;
    DoubleProperty m_orbitInclination;

    DoubleArrayProperty m_SAA_poly_lat;
    DoubleArrayProperty m_SAA_poly_lon;
    StringProperty m_xmlFiles;
    BooleanProperty m_aberrate;


};

// declare the service factories for the FluxSvc
//static SvcFactory<FluxSvc> a_factory;
//const ISvcFactory& FluxSvcFactory = a_factory;
DECLARE_SERVICE_FACTORY(FluxSvc);

static std::string default_source_library("$(FLUXXMLPATH)/source_library.xml");
static std::string default_dtd_file("$(FLUXXMLPATH)/source.dtd");

// ------------------------------------------------
// Implementation of the FluxSvc class
// ------------------------------------------------
/// Standard Constructor
FluxSvc::FluxSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc), m_currentFlux(0), m_insideSAA(false)
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
    declareProperty("SampleInterval", m_sampleInterval=1.0);
    declareProperty("OrbitInclination", m_orbitInclination=25.6);

    declareProperty("SAApolyLat"  , m_SAA_poly_lat);
    declareProperty("SAApolyLon"  , m_SAA_poly_lon);
    declareProperty("xmlListFile"    , m_xmlFiles="");
    declareProperty("EnableAberration", m_aberrate=false);

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

    // set orbit properties: needs to be done before EarthOrbit object gets created by GPS
    astro::EarthOrbit::set_inclination(m_orbitInclination);


    m_times.initialize(log);

    // set starting, or "launch" time for easy access by sources
    Spectrum::setStartTime(m_times.launch());


    status = serviceLocator()->queryInterface(IAppMgrUI::interfaceID(), (void**)&m_appMgrUI);

    // parse file with source library entries (consistent with obssim)
    if( !m_xmlFiles.value().empty() ) {

        std::string xmlFiles(m_xmlFiles.value() );
        facilities::Util::expandEnvVar(&xmlFiles);
        std::ifstream xmls(xmlFiles.c_str());
        if( !xmls.is_open() ){
            throw std::invalid_argument("File not found: " + xmlFiles);
        }
        while( ! xmls.eof()){
            std::string line; std::getline(xmls, line);
            if( line.empty() || line[0]=='#' ) continue; 
            m_source_lib.push_back(line);
        }
    }

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

    size_t nsaa =m_SAA_poly_lat.value().size();
    if( nsaa>0) {
        if( m_SAA_poly_lon.value().size() != nsaa ){
            log << MSG::ERROR <<"sizes of SAA arrays do not match"<< endreq;
            return StatusCode::FAILURE;
        }
        std::vector<std::pair<double,double> > saa_array;
        for( size_t i = 0; i< nsaa; ++i){
            saa_array.push_back( std::make_pair(m_SAA_poly_lat.value()[i], m_SAA_poly_lon.value()[i]) );
        }
        astro::EarthCoordinate::setSAAboundary( saa_array);
    }


    if( m_aberrate ){
        log << MSG::INFO << "Enabled generation of stellar aberration" << endreq;
    }
    astro::GPS::instance()->enableAberration(m_aberrate.value());


    //----------------------------------------------------------------
    // most of  the following cribbed from ToolSvc and ObjManager

    // look for a factory of an AlgTool that implements the IRegisterSource interface:
    // if found, make one and call the special method  

    // Manager of the AlgTool Objects
    /* IObjManager* objManager=0;             

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
*/

    // get a pointer to the tool service
    status = service( "ToolSvc", m_toolSvc, true );
    if (!status.isSuccess()) {
      log << MSG::ERROR << "Unable to get a handle to the tool service" << endmsg;
      return status;
    } else {
      log << MSG::DEBUG << "Got pointer to ToolSvc " << endmsg;
    }

  /* HMK Explicitly search for RegisterCRflux tool since objManager is
     no longer available in new Gaudi */

   IRegisterSource *ireg;
   status = m_toolSvc->retrieveTool("RegisterCRflux", ireg);
   if (status.isFailure())
       log << MSG::INFO << "No CRflux requested for this run" << endreq;
   else {
       log << MSG::INFO << "Found RegisterCRflux" << endreq;
       ireg->registerMe(this);
   }

   // In Interleave
   status = m_toolSvc->retrieveTool("RegisterSampledBackground", ireg);
   if (status.isSuccess()) {
       log << MSG::INFO << "Found RegisterSampledBackground" << endreq;
       ireg->registerMe(this);
   }

   // In userAlg
   status = m_toolSvc->retrieveTool("RegisterSource", ireg);
   if (status.isSuccess()) {
       log << MSG::INFO << "Found RegisterSource" << endreq;
       ireg->registerMe(this);
   }

    status = m_toolSvc->retrieveTool("FluxSvcRandom", m_randTool);
    if (status.isFailure()) 
        log << MSG::WARNING << "unable to create tool FluxSvcRandom" << endreq;

    m_fluxObs = new FluxObs();
    m_fluxObs->setFluxSvc(this);
    m_toolSvc->registerObserver(m_fluxObs);

    // attach an observer to be notified when orbital position changes
    // set callback to be notified when the position changes
    m_observer.setAdapter( new ActionAdapter<FluxSvc>
        (this, &FluxSvc::askGPS) );

    attachGpsObserver(&m_observer);

    astro::GPS::instance()->sampleintvl(m_sampleInterval.value()); // set sample interval for position updates
    return StatusCode::SUCCESS;
}
//------------------------------------------------------------------------
int FluxSvc::askGPS()
{
    astro::EarthCoordinate pos = GPS::instance()->earthpos();
    bool inside = pos.insideSAA();

    double curtime = GPS::instance()->time();

    if( m_insideSAA == inside) return 0; // no change

    MsgStream log( msgSvc(), name() );
    log << MSG::INFO;
    std::stringstream t;
        t<< (!inside? " leaving" : "entering") 
        << " SAA at "  << std::setprecision(10)<< curtime;
    log << t.str()    << endreq;
    m_insideSAA = inside;

    return 0; // can't be void in observer pattern
}


/// return pointer to the random engine
CLHEP::HepRandomEngine* FluxSvc::getRandomEngine(){
    return  CLHEP::HepRandom::getTheEngine();
};

// finalize
StatusCode FluxSvc::finalize ()
{
    StatusCode  status = StatusCode::SUCCESS;

    delete m_fluxMgr;
    return status;
}

/// Query interface
StatusCode FluxSvc::queryInterface(const InterfaceID& riid, void** ppvInterface)  {
    if ( IID_IFluxSvc.versionMatch(riid) )  {
        *ppvInterface = (IFluxSvc*)this;
    }else if (IRunable::interfaceID() == riid ) {
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

void FluxSvc::rootDisplay(std::vector<std::string> arguments){
    rootplot abc(arguments, m_fluxMgr);   
}


void FluxSvc::attachGpsObserver(Observer* anObserver)
{
    GPS::instance()->notification().attach( anObserver );
    GPS::instance()->notifyObservers(); // make sure everyone notified as observers are attached?
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


void FluxSvc::setPointingDirection(const astro::SkyDir& dir){
    astro::GPS::instance()->setPointingDirection(dir);
}


/// get the angular values of the satellite


CLHEP::HepRotation FluxSvc::transformToGlast(double time, GPS::CoordSystem index)const{
    return m_fluxMgr->transformToGlast(time,index);
}

CLHEP::HepRotation FluxSvc::transformGlastToGalactic(double time)const{
    return transformToGlast(time, GPS::CELESTIAL).inverse();
}

std::pair<double,double> FluxSvc::location(){
    return m_fluxMgr->location();
}
/// this sets the rocking mode in GPS.
std::vector<double> FluxSvc::setRockType(astro::GPS::RockType rockType, double rockAngle){
    return m_fluxMgr->setRockType(rockType, rockAngle);
}

std::vector<std::pair< std::string ,std::list<std::string> > > FluxSvc::sourceOriginList() const{
    return m_fluxMgr->sourceOriginList();
}

double FluxSvc::endruntime() {
    return m_times.end() ;
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
        IAlgManager::interfaceID(),
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

    GPS::instance()->notifyObservers(); // make sure all are in the

    // loop: will quit if either limit is set, and exceeded
    bool first(true);
    while( (m_evtMax==0  || m_evtMax>0 &&  eventNumber < m_evtMax)
        && (m_times.end()==0 || m_times.end()>0 && m_times.current() < m_times.end()) ) {

            double efrac =   (m_evtMax>0? 100.*eventNumber/m_evtMax: 0.0), 
                tfrac =   (m_times.end()>0? 100.*(m_times.current()-m_times.start())/(m_times.end()-m_times.start()) : 0.0) ;

            int percent_complete= static_cast<int>(  std::max( efrac, tfrac)  );
            if( percent_complete!=last_fraction){
                last_fraction=percent_complete;
                if( percent_complete<10 || percent_complete%10 ==0 || first){
                    first = false;

                /// time stamp for progress messages
		    facilities::Timestamp tstamp;
		    
                    log << MSG::INFO 
		        << " [" << tstamp.getString() << "]  "
                        //<<  std::setprecision(12)<< std::resetiosflags(4096) // scientific??
                        <<  std::setprecision(12)<< std::resetiosflags(std::ios::scientific) // scientific??
                        << percent_complete << "% complete: "
                        << " event "<< eventNumber<<",  time= " 
                        <<  m_times.current() << "= launch+ "
                        << (m_times.current()-m_times.launch()) << endreq;
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

    void FluxSvc::setAlignmentRotation(double qx, double qy, double qz, bool misalign)
    {
        m_fluxMgr->setAlignmentRotation(qx, qy, qz, misalign);
    }


    void FluxSvc::setFilterCone(std::vector<double> cone)
    {
        assert( cone.size()>2); // assume this already checked
        m_fluxMgr->setFilterCone(cone[0], cone[1], cone[2]);
    }


 
    void FluxSvc::setSAABoundary(const std::vector<std::pair<double, double> > & /*boundary*/)
    {
    }


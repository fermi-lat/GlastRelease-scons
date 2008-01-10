/** @file FluxAlg.cxx
@brief declaration and definition of the class FluxAlg

$Header$

*/

// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"
#include "GaudiKernel/SmartRefVector.h"

// Event for creating the McEvent stuff
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/Exposure.h"


// to write a Tree with entryosure info
#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include "facilities/Util.h"
#include "facilities/Observer.h"
#include "facilities/Timestamp.h"

#include "astro/GPS.h"

//flux
#include "FluxSvc/PointingInfo.h"
#include "FluxSvc/IFluxSvc.h"
#include "flux/IFlux.h"
#include "flux/Spectrum.h"

#include "flux/EventSource.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/Random.h"

#include <cassert>
#include <vector>
#include <iomanip>
#include <fstream>
#include <stdexcept>


// Include files
// Gaudi system includes
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"

using astro::GPS;

/** 
* \class FluxAlg
*
* \brief This is an Algorithm designed to get particle information 
* from FluxSvc and put it onto the TDS for later retrieval
* \author Toby Burnett
* 
* $Header$
*/

typedef HepGeom::Point3D<double>  HepPoint3D;
typedef HepGeom::Vector3D<double> HepVector3D;



class FluxAlg : public Algorithm {
public:
    FluxAlg(const std::string& name, ISvcLocator* pSvcLocator);
    double currentRate(){return m_currentRate;}

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();


private: 
    double m_currentRate;

    /// a single source
    StringProperty m_source_name;

    /// allow specification of a list of sources, which will be combined
    StringArrayProperty m_source_list;

    IFluxSvc*   m_fluxSvc;
    IFlux *     m_flux;

    unsigned int m_run;       // run number
    unsigned int m_sequence;  // sequence number

    IDataProviderSvc* m_eds;

    IParticlePropertySvc * m_partSvc;

    // the target area
    DoubleProperty m_area;
    DoubleProperty m_rocking_angle; // x-axis
    DoubleProperty m_rocking_angle_z; // z-axis

    std::map<int, int> m_counts; //! for measuring the total number generated per code.
    double m_initialTime;
    double m_currentTime;
    int m_SAAreject;

    PointingInfo m_pointing_info;

    
    StringArrayProperty m_pointingHistory;///< history file name and launch date

    StringProperty m_root_tree;
    BooleanProperty m_save_tuple; // set true to save
    BooleanProperty m_avoidSAA;


    INTupleWriterSvc* m_rootTupleSvc;

    ObserverAdapter< FluxAlg > m_observer; //obsever tag
    int askGPS(); // the function that will be called
    bool m_insideSAA; 

    IntegerProperty m_prescale;
    StringProperty m_source_info_filename;
    std::map<int, std::string> m_flux_names;
    void summary( std::ostream& log, std::string indent);
    DoubleArrayProperty m_misalignmentRotation;
    DoubleArrayProperty m_alignmentRotation;
    DoubleArrayProperty m_pointingDirection; ///< (ra, dec) for pointing
    DoubleProperty m_backoff; ///< backoff distance
    DoubleProperty m_zenithTheta; ///< set for zenith
    DoubleArrayProperty m_filterCone; ///< set parameters of a cone
    StringProperty(m_sourceListFile); ///< file name for list of sources (obssim compatible)
    BooleanProperty m_abort_on_exception; ///< what to do if an exception


};
//------------------------------------------------------------------------


static const AlgFactory<FluxAlg>  Factory;
const IAlgFactory& FluxAlgFactory = Factory;

//------------------------------------------------------------------------
//! ctor
FluxAlg::FluxAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator) , m_sequence(0), m_initialTime(0)
, m_SAAreject(0), m_insideSAA(false)
{
    
    // declare properties with setProperties calls
    declareProperty("source_name",  m_source_name="default");
    declareProperty("sources",     m_source_list);
    declareProperty("area",        m_area=6.0); // target area in m^2
    declareProperty("backoff",     m_backoff=2.0); //backoff distance in m

    declareProperty("rocking_angle", m_rocking_angle=0); // set non-zero to enable rocking

    declareProperty("PointingHistory",  m_pointingHistory); // doublet, filename and launch date

// deprecate these
    declareProperty("pointing_info_tree_name",  m_root_tree="");
    declareProperty("save_pointing_info",  m_save_tuple=false);

    declareProperty("AvoidSAA",   m_avoidSAA=false);
    declareProperty("Prescale",   m_prescale=1);
    declareProperty("source_info", m_source_info_filename="source_info.txt");
    declareProperty("misalignment", m_misalignmentRotation);
    declareProperty("alignment", m_alignmentRotation);
    declareProperty("pointingDirection", m_pointingDirection);
    declareProperty("zenithTheta", m_zenithTheta=-99);
    declareProperty("FilterCone",  m_filterCone);
    declareProperty("sourceListFile", m_sourceListFile="");
    declareProperty("abortOnException", m_abort_on_exception=false);

}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode FluxAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    GPS* gps = GPS::instance();

    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    // set target area for random point generation, and the backoff distance
    EventSource::totalArea(m_area);
    EventSource::s_backoff =m_backoff*1e3;// convert from m to mm

    // set pointing mode and associated rocking angle 
    if ( service("FluxSvc", m_fluxSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the FluxSvc!" << endreq;
        return StatusCode::FAILURE;
    }

    if( m_pointingDirection.value().size()==2 ) {
        // if direction set, ignore everything else!
        double ra(m_pointingDirection.value()[0]), dec(m_pointingDirection.value()[1]);
        m_fluxSvc->setPointingDirection( astro::SkyDir(ra, dec));
        log << MSG::INFO << "set to point at ra,dec= " << ra << ", "<<dec << endreq;

    }else if( m_zenithTheta>=0) {
        // if set from defalt, set to this value
        m_fluxSvc->setRockType(astro::GPS::EXPLICIT, m_zenithTheta);
        log << MSG::INFO << "set to zenith angle " << m_zenithTheta << " degrees" << endreq;

    }else {
        if(m_rocking_angle >0 ){
            //output to record the pointing settings
            //then this line sets the rocking type, as well as the rocking angle.
            m_fluxSvc->setRockType(GPS::ONEPERORBIT ,m_rocking_angle);
            log << MSG::INFO << "Once per orbit rocking Angle: " << m_rocking_angle << " degrees" << endreq;
        }
        
    
        //set the input file to be used as the pointing database, if used
        if(! m_pointingHistory.value().empty()){
            std::string filename(m_pointingHistory.value()[0]);
            facilities::Util::expandEnvVar(&filename);
            double offset = 0;
            bool horizontalflag(false);
            if( m_pointingHistory.value().size()>1){
                std::string field(m_pointingHistory.value()[1]);
                if(! field.empty() ) { // allow null string
                    facilities::Timestamp jt(m_pointingHistory.value()[1]);
                    offset = (astro::JulianDate(jt.getJulian())-astro::JulianDate::missionStart())*astro::JulianDate::secondsPerDay;
                }
            }

            if( m_pointingHistory.value().size()>2){
                std::string field(m_pointingHistory.value()[2]);
                horizontalflag =! field.empty();
            }
            log << MSG::INFO << "Loading Pointing History File : " << filename 
                << " with MET offset "<< offset <<  endreq;
            if( horizontalflag){
                log << MSG::INFO << "Will override x-direction to be horizontal"<<endreq;
            }

            gps->setPointingHistoryFile(filename, offset, horizontalflag);
        }
    }
    double current_time = gps->time(); // preserve time to protect against Pulsar, etc.

    // -------------- source name processing ------------------
    std::vector<std::string> sources(m_source_list.value()); // list from job options

    // parse file with source library entries (consistent with obssim)
    if( !m_sourceListFile.value().empty() ) {

        std::string fname(m_sourceListFile.value().c_str() );
        facilities::Util::expandEnvVar(&fname);
        std::ifstream file(fname.c_str());
        if( !file.is_open() ){
            throw std::invalid_argument("File not found: " + fname);
        }

        while( ! file.eof()){
            std::string line; std::getline(file, line);
            if( line.empty() || line[0]=='#' ) continue; 
            sources.push_back(line);
        }
    }

    // now, are there any in the list?

    if( !sources.empty()){
        log << MSG::INFO << "loading sources " << endreq;
        for(std::vector<std::string>::const_iterator it= sources.begin(); it!=sources.end(); ++it){
            log << MSG::INFO << "\t" << (*it) << endreq;
        }
        sc =  m_fluxSvc->compositeSource(sources, m_flux);
        if( sc.isFailure()) {
            log << MSG::ERROR << "Could not find one of the sources" << endreq;
            return sc;
        }


    }else{
        // no, try a single source for compatibility

        log << MSG::INFO << "loading source " << m_source_name << endreq;

        sc =  m_fluxSvc->source(m_source_name, m_flux);
        if( sc.isFailure()) {
            log << MSG::ERROR << "Could not find flux " << m_source_name << endreq;
            return sc;
        }
    }

    gps->time(current_time);  // restore time if it was modified
    gps->synch();             // and first notification of attached observers
    std::string title(m_flux->title()); if(title.length()>100) title = title.substr(0,100)+"...";
    log << MSG::INFO << "Source title: " << title << endreq;
    log << MSG::INFO << "        area: " << m_flux->targetArea() << endreq;
    log << MSG::INFO << "        rate: " << m_flux->rate() << endreq;
    if( m_prescale>1) {
        log << MSG::INFO << "    prescale: "<< m_prescale << endreq;
    }

    // check for (mis) alignment
    int alignment_pars ( m_alignmentRotation.value().size() );
    if (alignment_pars >0){

        double qx(m_alignmentRotation.value()[0]),qy(0),qz(0);
        if(alignment_pars >1) qy = m_alignmentRotation.value()[1];
        if(alignment_pars >2) qz = m_alignmentRotation.value()[2];
        m_fluxSvc->setAlignmentRotation(qx, qy, qz, false);
    }

    alignment_pars = ( m_misalignmentRotation.value().size() );
    if (alignment_pars >0){

        double qx(m_misalignmentRotation.value()[0]),qy(0),qz(0);
        if(alignment_pars >1) qy = m_misalignmentRotation.value()[1];
        if(alignment_pars >2) qz = m_misalignmentRotation.value()[2];
        m_fluxSvc->setAlignmentRotation(qx, qy, qz, true);
    }

    // check for filter cone
    if( m_filterCone.value().size()==3) {

        m_fluxSvc->setFilterCone(m_filterCone);
        
    }

    if ( service("ParticlePropertySvc", m_partSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the ParticlePropertySvc!" << endreq;
        return StatusCode::FAILURE;
    }

    // get a pointer to RootTupleSvc, use only if available 
    if( (service("RootTupleSvc", m_rootTupleSvc, true) ). isFailure() ) {
        log << MSG::WARNING << " RootTupleSvc is not available, will not write Pt tuple" << endreq;
        m_rootTupleSvc=0;
    }else if( !m_root_tree.value().empty() ) {
        
        m_pointing_info.setPtTuple(m_rootTupleSvc, m_root_tree.value());
    }


    // attach an observer to be notified when orbital position changes
        // set callback to be notified when the position changes
    m_observer.setAdapter( new ActionAdapter<FluxAlg>
        (this, &FluxAlg::askGPS) );

    m_fluxSvc->attachGpsObserver(&m_observer);
    gps->notifyObservers();
    return sc;
}
//------------------------------------------------------------------------
int FluxAlg::askGPS()
{
    astro::EarthCoordinate pos = GPS::instance()->earthpos();
    m_insideSAA = pos.insideSAA();

    return 0; // can't be void in observer pattern
}

//------------------------------------------------------------------------
//! process an event
StatusCode FluxAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    //
    // Purpose: have the flux service create parameters of an incoming particle 
    // if nothing has changed, then use the existing m_flux,
    // but if the "current" IFlux is not the same as the one we have now,
    // then change our m_flux pointer to be the new one.
    // Output:  a staturCode to ensure the function executed properly.
    m_flux = m_fluxSvc->currentFlux();

    // check the current random number seed
    int seed = CLHEP::HepRandom::getTheSeed();
    log << MSG::DEBUG << "random seed: " << seed << endreq;

    std::string particleName;
    if( m_avoidSAA) m_SAAreject--;
    int count = m_prescale;
    do{ // loop if we are rejecting particles generated during SAA
        // also do prescale here
        try {
            bool valid =m_flux->generate();
            if( !valid) {
                log << MSG::ERROR << "Ran out of valid sources, aborting" << endreq;
                return StatusCode::FAILURE;
            }
        } catch( const std::exception & e) {
            std::cout << "FluxAlg caught exception, aborting this event " << e.what() 
                << "\n\tprocessing source " << m_flux->particleName() 
                << "\n\tcurrent time: " <<GPS::instance()->time() << std::endl;
            if( m_abort_on_exception ){
                return StatusCode::FAILURE;
            }
            setFilterPassed( false); // should go to clocks.
            return StatusCode::SUCCESS;
            
        }
        particleName = m_flux->particleName();

        //if it's a clock then ExposureAlg will take care of it, and no othe algorithms should care about it.
        if(particleName == "TimeTick" || particleName == "Clock"){
            m_pointing_info.set();

            setFilterPassed( false );
            return sc;
        }
        if(m_insideSAA && m_avoidSAA.value() ){
            double time(GPS::instance()->time()), 
                endtime( m_fluxSvc->endruntime() );
            if( time >endtime ){
                log << MSG::INFO << "Ran out of time while in SAA"<< endreq;
                setFilterPassed( false );
                break;  //return sc;
            }
                
            ++m_SAAreject;
        }
        else break;
    } while(m_insideSAA && m_avoidSAA.value() || --count>0);

    Hep3Vector p = m_flux->launchPoint();
    Hep3Vector d = m_flux->launchDir();

    double ke = m_flux->energy(); // kinetic energy in MeV

    //here's where we get the particleID and mass for later.
    // Note that the Gaudi particle table now only has p+ and n0 for proton and neutron: 
    if( particleName=="p" || particleName=="proton") particleName="p+";
    if( particleName=="neutron") particleName="n0";
    ParticleProperty* prop = m_partSvc->find(particleName);

//    if( prop==0 && particleName=="He" ){
//        // If He didn't work (mystery!) try alpha instead
//        prop = m_partSvc->find("alpha");
//    }
    if( prop==0) {
        log << MSG::ERROR << "Particle name " << particleName << " not found by particle properties" << endreq;
        return StatusCode::FAILURE;
    }

    int partID = prop->jetsetID(); // same as stdhep id

    log << MSG::DEBUG ;
    if( log.isActive()){
        log<< particleName << ", flux("<<m_flux->name() << ") "
        << "(" << m_flux->energy()
        << " MeV), Launch: " 
        << "(" << p.x() <<", "<< p.y() <<", "<<p.z()<<")" 
        << " mm, Dir " 
        << "(" << d.x() <<", "<< d.y() <<", "<<d.z()<<")" ;
    }
    log << endreq;


    // Here the TDS is prepared to receive hits vectors
    // Check for the MC branch - it will be created if it is not available
    Event::MCEvent* mch = 0;

    SmartDataPtr<Event::MCEvent> mcheader(eventSvc(), EventModel::MC::Event);
    if (mcheader == 0) {
        sc=eventSvc()->registerObject(EventModel::MC::Event , mch= new Event::MCEvent);
        mch->initialize(0,0,m_sequence, m_flux->time(), m_flux->name());
        if(sc.isFailure()) {
            log << MSG::WARNING << EventModel::MC::Event  <<" could not be registered on data store" << endreq;
            delete mch;
            return sc;
        }

    }else {
        mch = mcheader;
    }

    mch->initialize(mch->getRunNumber(), m_flux->numSource(), mch->getSequence(), m_flux->time());
    mch->setSourceName(m_flux->name());

    Event::McParticleCol* pcol = new Event::McParticleCol;
    sc = eventSvc()->registerObject(EventModel::MC::McParticleCol, pcol);
    if( sc.isFailure()) {

        log << MSG::ERROR << "Could not Register "<< EventModel::MC::McParticleCol << endreq;

        return sc;
    }
    Event::McParticle * parent= new Event::McParticle;
    pcol->push_back(parent);

    double mass = prop->mass() , 
        energy = (ke+mass),
        momentum=sqrt(energy*energy - mass*mass); 
    CLHEP::HepLorentzVector pin(d*momentum,energy);

    // This parent particle decay at the start in the first particle, 
    // so initial momentum and final one are the same
    parent->initialize(parent, partID, 
        Event::McParticle::PRIMARY,
        pin,p, m_flux->name());
    parent->finalize(pin, p);

    // get the event header to set the time
    Event::EventHeader* h = 0; 

    SmartDataPtr<Event::EventHeader> header(eventSvc(), EventModel::EventHeader);
    if(0==header) {
        // not already there: try to register instead
        //sc = eventSvc()->registerObject(EventModel::EventHeader, h=new Event::EventHeader);
        sc = eventSvc()->registerObject(EventModel::EventHeader, EventModel::EventHeader, h=new Event::EventHeader);
        if( sc.isFailure()) {
            log << MSG::WARNING << " could not find or register the event header" << endreq;
        }
    } else{  h = header;
    }

    m_currentTime=m_flux->time();

    // is this proper here?
    m_pointing_info.set();
    
    // put pointing stuff into the root tree
    if( m_rootTupleSvc!=0 && !m_root_tree.value().empty()){
        m_rootTupleSvc->storeRowFlag(this->m_root_tree.value(), m_save_tuple);
    }

    if( m_initialTime==0) m_initialTime=m_currentTime;
    h->setTime(m_currentTime);
    m_run = h->run();  // save
    int numEvents = ++m_sequence;
    m_counts[m_flux->numSource()]++; // update count
    m_flux_names[m_flux->numSource()]= m_flux->name(); // save (or resave!) name 

    m_currentRate=numEvents/(m_currentTime-m_initialTime);
    return StatusCode::SUCCESS;
}

namespace { 
    // define some static variables for rootupleSVc to get after we are deleted
    unsigned int run, sequence;
    double initialTime, currentTime;
}
//------------------------------------------------------------------------
//! clean up, summarize
StatusCode FluxAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    static bool done = false;
    if( done  ) return sc;
    done=true;

    if( m_rootTupleSvc!=0 ){
        // create the jobinfo tuple: copy to statics
        run = m_run;
        sequence=m_sequence + (m_SAAreject>0? m_SAAreject: 0 );
        initialTime=m_initialTime;
        currentTime=m_currentTime;
        m_rootTupleSvc->addItem("jobinfo", "run", &run);
        m_rootTupleSvc->addItem("jobinfo", "generated", &sequence);
        m_rootTupleSvc->addItem("jobinfo", "start", &initialTime);
        m_rootTupleSvc->addItem("jobinfo", "stop",  &currentTime);
    }
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "Computed Rate: "<< currentRate() << " Hz" ;
    summary(log.stream(), "\n\t\t\t");
    log  << endreq;

    if( !m_source_info_filename.value().empty() ){
        std::ofstream infofile(m_source_info_filename.value().c_str());
        summary( infofile, std::string("\n") );
    }

    
    if( m_avoidSAA && m_SAAreject>0 ){
        log << "\t\tRejected by SAA: " << m_SAAreject << endreq;
            log << "\t\t(note that this may invalidate the rate calculation)" << endreq;
    }
    return sc;
}

void FluxAlg::summary( std::ostream& log, std::string indent)
{
    log << indent << " Source ID   Source Name                  counts";
    for(std::map<int,int>::const_iterator im=m_counts.begin(); im !=m_counts.end(); ++im) {
        log << indent
            << std::setw(10) <<im->first 
            << "   "  << std::setw(25) << std::left<< m_flux_names[im->first]
            << std::setw(10)<< std::right << im->second;
    }
    log << std::endl;
}


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

//flux
#include "FluxSvc/PointingInfo.h"
#include "FluxSvc/IFluxSvc.h"
#include "flux/IFlux.h"
#include "flux/Spectrum.h"
#include "flux/SpectrumFactory.h"

#include "flux/EventSource.h"
#include "flux/CompositeSource.h"

#include "CLHEP/Vector/LorentzVector.h"

#include <cassert>
#include <vector>


// Include files
// Gaudi system includes
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"


/** 
* \class FluxAlg
*
* \brief This is an Algorithm designed to get particle information 
* from FluxSvc and put it onto the TDS for later retrieval
* \author Toby Burnett
* 
* $Header$
*/


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


    UnsignedIntegerProperty m_run;      // run number
    unsigned int m_sequence;  // sequence number
    StringProperty m_pointing_history_input_file;

    IDataProviderSvc* m_eds;

    IParticlePropertySvc * m_partSvc;

    // the target area
    DoubleProperty m_area;
    IntegerProperty m_pointing_mode;
    DoubleProperty m_rocking_angle; // x-axis
    DoubleProperty m_rocking_angle_z; // z-axis

    std::map<int, int> m_counts; //! for measuring the total number generated per code.
    TimeStamp m_initialTime;
    int m_SAAreject;

    PointingInfo m_pointing_info;
    StringProperty m_root_tree;
    BooleanProperty m_save_tuple; // set true to save
    BooleanProperty m_avoidSAA;


    INTupleWriterSvc* m_rootTupleSvc;

    ObserverAdapter< FluxAlg > m_observer; //obsever tag
    int askGPS(); // the function that will be called
    bool m_insideSAA; 

    IntegerProperty m_prescale;

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
    declareProperty("MCrun",        m_run=100);
    declareProperty("area",        m_area=6.0); // target area in m^2
    declareProperty("pointing_mode", m_pointing_mode=0);
    declareProperty("rocking_angle", m_rocking_angle=0); // in degrees
    declareProperty("rocking_angle_z", m_rocking_angle_z=0); // in degrees
    declareProperty("pointing_history_input_file",  m_pointing_history_input_file="");

    declareProperty("pointing_info_tree_name",  m_root_tree="MeritTuple");
    declareProperty("save_pointing_info",  m_save_tuple=false);
    declareProperty("AvoidSAA",   m_avoidSAA=false);
    declareProperty("Prescale",   m_prescale=1);


}
//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode FluxAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;

    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    // set target area for random point generation
    EventSource::totalArea(m_area);

    // set pointing mode and associated rocking angle 
    if ( service("FluxSvc", m_fluxSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the FluxSvc!" << endreq;
        return StatusCode::FAILURE;
    }

    //this line sets the explicit rocking angles to be used IF the 
    //rocking type is "explicit." Note that it uses the m_rocking_angle also.
    m_fluxSvc->setExplicitRockingAngles(m_rocking_angle*M_PI/180,m_rocking_angle_z*M_PI/180);
    //m_fluxSvc->setRockingAngle(m_rocking_angle);

    //then this line sets the rocking type, as well as the rocking angle.
    m_fluxSvc->setRockType(m_pointing_mode,m_rocking_angle);


    //set the input file to be used as the pointing database, if used
    if(! m_pointing_history_input_file.value().empty() ){
        m_fluxSvc->setPointingHistoryFile(m_pointing_history_input_file.value().c_str());
    }

    //now make the screen outputs showing what modes were called:
    if(m_pointing_mode){
        //output to record the pointing settings
        log << MSG::INFO << "rocking Mode: " << m_pointing_mode << endreq;
        log << MSG::INFO << "rocking Angle: " << m_rocking_angle << " degrees" << endreq;
    }
    if(! m_pointing_history_input_file.value().empty()){
        log << MSG::INFO << "Pointing History File used, file name: " << m_pointing_history_input_file.value() << endreq;
        if(m_pointing_mode){
            //problem - two kinds of pointing are being used!
            log << MSG::WARNING << "Pointing History and rocking mode both specified!" << endreq;
        }else{
            //no rocking called, mode=0.  say so.
            log << MSG::INFO << "no rocking mode called, none required." << endreq;
        }
    }

    //finally, if no rocking or pointing of any kind has been set:
    if(!m_pointing_mode && m_pointing_history_input_file.value().empty()){
        log << MSG::INFO << "no rocking mode called, rocking disabled for this run." << endreq;
    }

    if( !m_source_list.value().empty()){
        log << MSG::INFO << "loading sources " << endreq;
        std::vector<std::string> sources=m_source_list.value();
        for(std::vector<std::string>::const_iterator it= sources.begin(); it!=sources.end(); ++it){
            log << MSG::INFO << "\t" << (*it) << endreq;
        }
        sc =  m_fluxSvc->compositeSource(sources, m_flux);
        if( sc.isFailure()) {
            log << MSG::ERROR << "Could not find one of the sources" << endreq;
            return sc;
        }


    }else{
        log << MSG::INFO << "loading source " << m_source_name << endreq;

        sc =  m_fluxSvc->source(m_source_name, m_flux);
        if( sc.isFailure()) {
            log << MSG::ERROR << "Could not find flux " << m_source_name << endreq;
            return sc;
        }
    }
    std::string title(m_flux->title()); if(title.length()>100) title = title.substr(0,100)+"...";
    log << MSG::INFO << "Source title: " << title << endreq;
    log << MSG::INFO << "        area: " << m_flux->targetArea() << endreq;
    log << MSG::INFO << "        rate: " << m_flux->rate() << endreq;
    if( m_prescale>1) {
        log << MSG::INFO << "    prescale: "<< m_prescale << endreq;
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

    std::string particleName;
    m_SAAreject--;
    int count = m_prescale;
    do{ // loop if we are rejecting particles generated during SAA
        // also do prescale here
        m_flux->generate();
        particleName = m_flux->particleName();

        //if it's a clock then ExposureAlg will take care of it, and no othe algorithms should care about it.
        if(particleName == "TimeTick" || particleName == "Clock"){
            setFilterPassed( false );
            return sc;
        }
        ++m_SAAreject;
    } while(m_insideSAA && m_avoidSAA.value() || --count>0);

    HepPoint3D p = m_flux->launchPoint();
    HepVector3D d = m_flux->launchDir();
    double ke = m_flux->energy(); // kinetic energy in MeV

    //here's where we get the particleID and mass for later.
    if( particleName=="p") particleName="proton";
    ParticleProperty* prop = m_partSvc->find(particleName);

    if( prop==0) {
        log << MSG::ERROR << "Particle name " << particleName << " not found by particle properties" << endreq;
        return StatusCode::FAILURE;
    }

    int partID = prop->jetsetID(); // same as stdhep id

    log << MSG::DEBUG ;
    if( log.isActive()){
        log<< particleName
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
        mch->initialize(0,0,m_sequence, m_flux->time());
        if(sc.isFailure()) {
            log << MSG::WARNING << EventModel::MC::Event  <<" could not be registered on data store" << endreq;
            delete mch;
            return sc;
        }

    }else mch = mcheader;


    mch->initialize(mch->getRunNumber(), m_flux->numSource(), mch->getSequence(), m_flux->time());



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
    HepLorentzVector pin(d*momentum,energy);

    // This parent particle decay at the start in the first particle, 
    // so initial momentum and final one are the same
    parent->initialize(parent, partID, 
        Event::McParticle::PRIMARY,
        pin,p);
    parent->finalize(pin, p);

    // get the event header to set the time
    Event::EventHeader* h = 0; 

    SmartDataPtr<Event::EventHeader> header(eventSvc(), EventModel::EventHeader);
    if(0==header) {
        // not already there: try to register instead
        sc = eventSvc()->registerObject(EventModel::EventHeader, h=new Event::EventHeader);
        if( sc.isFailure()) {
            log << MSG::WARNING << " could not find or register the event header" << endreq;
        }
    } else{  h = header;
    }

    TimeStamp currentTime=m_flux->time();

    m_pointing_info.set(currentTime);

    // Here the TDS receives the exposure data
    Event::ExposureCol* exposureDBase = new Event::ExposureCol;
    sc=eventSvc()->registerObject(EventModel::MC::ExposureCol , exposureDBase);
    if(sc.isFailure()) {
        log << MSG::ERROR << EventModel::MC::ExposureCol  
            <<" could not be entered into existing data store" << endreq;
        return sc;
    }
    exposureDBase->push_back(m_pointing_info.forTDS());

    
    // put pointing stuff into the root tree
    if( m_rootTupleSvc!=0 && !m_root_tree.value().empty()){
        m_rootTupleSvc->storeRowFlag(this->m_root_tree.value(), m_save_tuple);
    }

    if( m_initialTime==0) m_initialTime=currentTime;
    h->setTime(currentTime);
    int numEvents = ++m_sequence;
    m_counts[m_flux->numSource()]++; // update count

    m_currentRate=numEvents/(currentTime-m_initialTime);
    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode FluxAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    static bool done = false;
    if( done || m_counts.empty() ) return sc;
    done=true;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "Computed Rate: "<< currentRate() << " Hz"
        << "\n\t\t\t  Source ID    counts";
    for(std::map<int,int>::const_iterator im=m_counts.begin(); im !=m_counts.end(); ++im) {
        log << "\n\t\t\t\t" << im->first << "\t" << im->second;
    }
    log  << endreq;

    if( m_avoidSAA){
        log << "\t\tRejected by SAA: " << m_SAAreject << endreq;
            log << "\t\t(note that this may invalidate the rate calculation)" << endreq;
    }
    return sc;
}


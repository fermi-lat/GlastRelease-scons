/** 
* @file ExposureAlg.cxx
* @brief Definition and implementation of class ExposureAlg
*
*  $Header$
*/

// Include files
// Gaudi system includes
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"

#include "astro/SkyDir.h"
#include "astro/EarthCoordinate.h"
#include "astro/SolarSystem.h"
#include "astro/EarthOrbit.h"
// Event for creating the McEvent stuff
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/Exposure.h"

//flux
#include "FluxSvc/IFluxSvc.h"
#include "FluxSvc/IFlux.h"
#include "astro/GPS.h"

// to write a Tree with entryosure info
#include "ntupleWriterSvc/INTupleWriterSvc.h"
#include "facilities/Util.h"

#include <cassert>
#include <vector>
#include <fstream>

/** 
* \class ExposureAlg
*
* \brief This is an Algorithm designed to get information about LAT
* position,  from FluxSvc and use it to put information onto the TDS about
* LAT pointing and location characteristics, effectively generating the D2 database.  The "TimeTick" 
* Spectrum is included (and can be used in jobOptions with this algorithm) in order to provide a constant time reference.
*
* \author Sean Robinson
* 
* $Header$
*/
class ExposureAlg : public Algorithm {
public:
    ExposureAlg(const std::string& name, ISvcLocator* pSvcLocator);

    //stuff that an Algorithm needs.
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

        /** nested class to manage extra TTrees 
    */
    class TTree { 
    public:
        /**
        */
        TTree(     INTupleWriterSvc* rootTupleSvc, 
            std::string treeName,  
            std::vector<const char* >leaf_names)
            :m_name(treeName), m_leafNames(leaf_names)
        {
            m_values.resize(leaf_names.size());
            int i=0;
            for( std::vector<const char*>::const_iterator  it = leaf_names.begin();
                it!=leaf_names.end(); 
                ++it){
                    rootTupleSvc->addItem(treeName,  *it,  &m_values[i++]);
                }
        }
        void fill(int n, double value){ m_values[n]=value;};
        /** called from gui printer */
        void printOn(std::ostream& out)const {
            out << "Tree "<< m_name << std::endl;
            std::vector<double>::const_iterator dit=m_values.begin();
            for( std::vector<const char*> ::const_iterator nit=m_leafNames.begin(); nit!=m_leafNames.end(); ++nit,++dit){
                out << std::setw(15) << *nit << "  " <<  *dit << std::endl;
            }
        }
        /// data values: name of tree, leaves, and the values
        std::string m_name;
        std::vector<const char*> m_leafNames;
        std::vector< double> m_values;
    };
private: 
    TTree* m_pointing_tree;;


    double m_lasttime; //time value to hold time between events;
    StringProperty m_source_name;
    StringProperty m_pointing_history_output_file;
    StringProperty m_pointing_history_input_file;
    StringProperty m_root_tree;

    IFluxSvc*   m_fluxSvc;
    IFlux *     m_flux;

    std::ostream* m_out;  //for output that looks like the stuff from the astro orbit model test.
    int         m_tickCount; // number of ticks processed
    INTupleWriterSvc* m_rootTupleSvc;;


};
//------------------------------------------------------------------------

static const AlgFactory<ExposureAlg>  Factory;
const IAlgFactory& ExposureAlgFactory = Factory;

//------------------------------------------------------------------------
//! ctor
ExposureAlg::ExposureAlg(const std::string& name, ISvcLocator* pSvcLocator)
: Algorithm(name, pSvcLocator)
, m_pointing_tree(0)
, m_out(0)
, m_lasttime(0)
, m_tickCount(0)
{
    // declare properties with setProperties calls
    declareProperty("source_name",  m_source_name="default");
    declareProperty("pointing_history_output_file",  m_pointing_history_output_file="");
    declareProperty("pointing_history_input_file",  m_pointing_history_input_file="");
    declareProperty("root_tree",  m_root_tree="");

}

//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode ExposureAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;

    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    if ( service("FluxSvc", m_fluxSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the FluxSvc!" << endreq;
        return StatusCode::FAILURE;
    }

    // access the properties of FluxSvc
    IProperty* propMgr=0;
    sc = serviceLocator()->service("FluxSvc", propMgr );
    if( sc.isFailure()) {
        log << MSG::ERROR << "Unable to locate PropertyManager Service" << endreq;
        return sc;
    }

    DoubleProperty startTime("StartTime",0);
    sc = propMgr->getProperty( &startTime );
    if (sc.isFailure()) return sc;

    m_lasttime = startTime.value();

    //set the input file to be used as the pointing database
    if(! m_pointing_history_input_file.value().empty() ){
        std::string fileName(m_pointing_history_input_file.value());
        facilities::Util::expandEnvVar(&fileName);
        m_fluxSvc->setPointingHistoryFile(fileName.c_str());
    }

    //set the output file (pointing information) to be written.
    if(! m_pointing_history_output_file.value().empty() ){
        std::string fileName(m_pointing_history_output_file.value());
        facilities::Util::expandEnvVar(&fileName);
        m_out = new std::ofstream(fileName.c_str());
    }
    if(!m_root_tree.value().empty() ){
        // get a pointer to RootTupleSvc 
        if( (sc = service("RootTupleSvc", m_rootTupleSvc, true) ). isFailure() ) {
            log << MSG::ERROR << " failed to get the RootTupleSvc" << endreq;
            return sc;
        }
        // now define the root tuple
        std::vector<const char* > names;
        const char * point_info_name[] = {"time","lat","lon","alt","posx","posy","posz","rax","decx","raz","decz"};
        for( int i = 0; i< (int)(sizeof(point_info_name)/sizeof(void*)); ++i){ 
            names.push_back(point_info_name[i]); }

        m_pointing_tree = new TTree( m_rootTupleSvc,  std::string(m_root_tree),  names);

    }

    return sc;
}


//------------------------------------------------------------------------
//! process an event
StatusCode ExposureAlg::execute()
{
    using namespace astro;
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    double currentTime;
    //-----------------------------------------------------------------------

    Event::McParticleCol* pcol = 0;
    eventSvc()->retrieveObject("/Event/MC/McParticleCol",(DataObject *&)pcol);
    //only make a new source if one does not already exist.
    if(pcol==0){
        //FluxAlg didn't do anything.  proceed.

        //if the flux had changed, something changed the source type.
        if(m_fluxSvc->currentFlux() == m_flux){
            m_flux->generate();
            currentTime = m_flux->time();
        }else{
            m_flux = m_fluxSvc->currentFlux();
            m_flux->generate();
            currentTime = m_flux->time();
        }
    }else{
        //FluxAlg is taking care of the particle, so do nothing but get the current time.
        IFlux * f = m_fluxSvc->currentFlux();
        if( f!=0)   currentTime =f->time();
        else{
            // get time from the McHeader 
                SmartDataPtr<Event::MCEvent>     mcheader(eventSvc(),    EventModel::MC::Event);
               currentTime= mcheader->time();
        }
    }   
    // is this a special "particle" designed to save the information?
    bool tick =    m_fluxSvc->currentFlux()!=0 && m_fluxSvc->currentFlux()->particleName() == "TimeTick" ;

    //by now, we should know that we have the appropriate particle to make a Exposure with.
    // here we get the time characteristics



    //NOTE: this gets an interval from the last time that a TimeTick particle came to this one.
    //in other words, the timeTick particles define the beginning and ends of intervals.
    double intrvalstart = m_lasttime;
    //then get the end of the interval (this can be made into whatever).
    double intrvalend = currentTime;
    //..and reset the time of this event to be the "last time" for next time.
    if(tick) m_lasttime = intrvalend;


    // The GPS singleton has current time and orientation
    GPS* gps = GPS::instance();

    // evaluating pointing stuff for previous tick, or current time if not a tick
    double time = tick? intrvalstart : currentTime;
    gps->getPointingCharacteristics(time); // sets time for other functions
    Hep3Vector location = gps->position(time); // special, needs its own time

    // hold onto the cartesian location of the LAT
    double posx = location.x(), 
        posy = location.y(), 
        posz = location.z(); 

    // directions of the x and z axes, and the zenith
    double 
        rax =   gps->RAX(),        decx =  gps->DECX(),
        raz =   gps->RAZ(),        decz =  gps->DECZ(),
        razenith = gps->RAZenith(),deczenith = gps->DECZenith();
    //uncomment for debug check double check=astro::SkyDir(rax, decx)().dot(astro::SkyDir(raz, decz)());

    EarthOrbit orb; //for the following line - this should have a better implementation.
    double julianDate = orb.dateFromSeconds(m_lasttime);
    EarthCoordinate earthpos(location,julianDate);
    double lat = earthpos.latitude();
    double lon = earthpos.longitude();
    double alt = earthpos.altitude();

#if 0 // example code, not needed if we don't put sun and moon, SAA info in yet (which are redundant)
    bool SAA = earthpos.insideSAA();
    SolarSystem sstm;

    double ramoon =  sstm.direction(astro::SolarSystem::Moon,julianDate).ra();
    double decmoon = sstm.direction(astro::SolarSystem::Moon,julianDate).dec();
    double rasun =   sstm.direction(astro::SolarSystem::Sun,julianDate).ra();
    double decsun =  sstm.direction(astro::SolarSystem::Sun,julianDate).dec();
    SkyDir sunDir(rasun,decsun);
    //Rotation galtoglast(m_fluxSvc->transformGlastToGalactic(currentTime).inverse);
    sunDir()=(m_fluxSvc->transformGlastToGalactic(intrvalstart).inverse())*sunDir();

#endif

    // Here the TDS receives the exposure data
    Event::ExposureCol* exposureDBase = new Event::ExposureCol;
    sc=eventSvc()->registerObject(EventModel::MC::ExposureCol , exposureDBase);
    if(sc.isFailure()) {
        log << MSG::ERROR << EventModel::MC::ExposureCol  
            <<" could not be entered into existing data store" << endreq;
        return sc;
    }

    Event::Exposure* entry = new Event::Exposure;
    exposureDBase->push_back(entry);
    entry->init(tick?intrvalstart:currentTime,lat,lon,alt,posx,posy,posz,rax,decx,raz,decz);

    //now, only do the rest of this algorithm if we have a timetick particle.

    if(    !tick ) return StatusCode::SUCCESS;

    log << MSG::DEBUG ;
    if(log.isActive()){

        // now we'll retreive the data from the TDS as a check.
        Event::ExposureCol* elist = 0;
        eventSvc()->retrieveObject("/Event/MC/ExposureCol",(DataObject *&)elist);

        Event::ExposureCol::iterator curEntry = (*elist).begin();
        //some test output - to show that the data got onto the TDS
        (*curEntry)->fillStream(log.stream());
    }
    log << endreq;
#if 0 // another example
    SkyDir curDir(raz,decz);
    SkyDir xDir(rax,decx);
#endif

    //and here's the file output.
    if( m_out !=0) {
        std::ostream& out = *m_out;
        out <<  std::setprecision(14) 
            <<intrvalstart <<'\t';
        out << std::setw(10) << std::setprecision(10)
            << posx<<'\t';
        out<< posy<<'\t';
        out<< posz<<'\t';
        out<<raz<<'\t';
        out<<decz<<'\t';
        out<<rax<<'\t';
        out<<decx<<'\t';
        out<<razenith<<"\t";
        out<<deczenith <<'\t';
        out<<lon <<'\t';
        out<<lat <<'\t';
        out<<alt << std::endl;

    }
    //-----------------------------------------------
    if( m_pointing_tree !=0 ) {
    int n= 0;
        m_pointing_tree->fill(n++,entry->intrvalstart());
        m_pointing_tree->fill(n++,entry->lat());
        m_pointing_tree->fill(n++,entry->lon());
        m_pointing_tree->fill(n++,entry->alt());
        m_pointing_tree->fill(n++,entry->posX());
        m_pointing_tree->fill(n++,entry->posY());
        m_pointing_tree->fill(n++,entry->posZ());
        m_pointing_tree->fill(n++,entry->RAX());
        m_pointing_tree->fill(n++,entry->DECX());
        m_pointing_tree->fill(n++,entry->RAZ());
        m_pointing_tree->fill(n++,entry->DECZ());
        m_rootTupleSvc->storeRowFlag(true);
    }
    //----------------------------------------------
    
    m_tickCount++;
    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode ExposureAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "Processed " << m_tickCount << " ticks" << endreq;
    delete m_out;

    return sc;
}


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

//flux
#include "FluxSvc/IFluxSvc.h"
#include "flux/IFlux.h"
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


        // the root tuple, a la FT2
        // http://glast.gsfc.nasa.gov/ssc/dev/fits_def/definitionFT2.html

    class FT2entry {
    public:
        FT2entry(INTupleWriterSvc* tuple, const std::string& tname);

        double start, stop;
        float sc_position[3];
        float lat_geo, lon_geo;
        float rad_geo;
        float ra_zenith, dec_zenith;
        float ra_scz, dec_scz;
        float ra_scx, dec_scx;
        float livetime;

        void begin(double start_time);
        void finish(double stop_time, double live);
    };
private: 

    double m_lasttime; //time value to hold time between events;
    double m_initial_time; 

    StringProperty m_root_tree;
    StringProperty m_pointing_history_input_file;

    IFluxSvc*   m_fluxSvc;

    int         m_tickCount; // number of ticks processed
    INTupleWriterSvc* m_rootTupleSvc;;

    FT2entry* m_history;
};
//------------------------------------------------------------------------

static const AlgFactory<ExposureAlg>  Factory;
const IAlgFactory& ExposureAlgFactory = Factory;

void ExposureAlg::FT2entry::begin(double time)
{
    // Purpose: save the current status until the next tick
    using namespace astro;
 
    start = time;
    // The GPS singleton has current time and orientation
    GPS* gps = GPS::instance();
    gps->getPointingCharacteristics(time); // sets time for other functions
    Hep3Vector location = 1.e3* gps->position(time); // special, needs its own time
    
    // cartesian location of the LAT
    sc_position[0] = location.x();
    sc_position[1] = location.y();
    sc_position[2] = location.z(); 

    ra_zenith = gps->RAZenith();
    dec_zenith= gps->DECZenith();
    ra_scx =    gps->RAX();
    dec_scx =   gps->DECX();

    ra_scz =    gps->RAZ();
    dec_scz =   gps->DECZ();
    //uncomment for debug check double check=astro::SkyDir(rax, decx)().dot(astro::SkyDir(raz, decz)());

    lat_geo = gps->lat(); 
    lon_geo = gps->lon(); 
    rad_geo = gps->altitude(); 
}
void ExposureAlg::FT2entry::finish(double stop_time, double live)
{
    stop = stop_time;
    livetime = live;
}
ExposureAlg::FT2entry::FT2entry(INTupleWriterSvc* tuple, const std::string& tname)
{
    if( tuple==0 ) return;
    tuple->addItem(tname, "start",  &start);
    tuple->addItem(tname, "stop",   &stop);
    tuple->addItem(tname, "sc_position[3]", sc_position);
    tuple->addItem(tname, "lat_geo", &lat_geo);
    tuple->addItem(tname, "lon_geo", &lon_geo);
    tuple->addItem(tname, "rad_geo", &rad_geo);
    tuple->addItem(tname, "ra_zenith", &ra_zenith);
    tuple->addItem(tname, "dec_zenith",&dec_zenith);
    tuple->addItem(tname, "ra_scz", &ra_scz);
    tuple->addItem(tname, "dec_scz", &dec_scz);
    tuple->addItem(tname, "ra_scx",   &ra_scx);
    tuple->addItem(tname, "dec_scx", &dec_scx);
    tuple->addItem(tname, "livetime", &livetime);
}
//------------------------------------------------------------------------
//! ctor
ExposureAlg::ExposureAlg(const std::string& name, ISvcLocator* pSvcLocator)
: Algorithm(name, pSvcLocator)
, m_lasttime(0)
, m_tickCount(0)
{
    // declare properties with setProperties calls
    declareProperty("root_tree",  m_root_tree="pointing_history");
    declareProperty("pointing_history_input_file",  m_pointing_history_input_file="");

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

    //set the input file to be used as the pointing database
    if(! m_pointing_history_input_file.value().empty() ){
        std::string fileName(m_pointing_history_input_file.value());
        facilities::Util::expandEnvVar(&fileName);
        GPS::instance()->setPointingHistoryFile(fileName);
    }


    DoubleProperty startTime("StartTime",0);
    sc = propMgr->getProperty( &startTime );
    if (sc.isFailure()) return sc;

    m_initial_time =m_lasttime = startTime.value();

    // get a pointer to RootTupleSvc 
    if( (sc = service("RootTupleSvc", m_rootTupleSvc, true) ). isFailure() ) {
            log << MSG::ERROR << " failed to get the RootTupleSvc" << endreq;
            return sc;
    }

    m_history = new FT2entry(m_rootTupleSvc, m_root_tree.value());

    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode ExposureAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    m_lasttime = m_fluxSvc->currentFlux()->time();

    log << MSG::DEBUG ;
    if( log.isActive() ){
        std::stringstream t;
        t << "tick at " << std::setprecision(10)
        << m_lasttime - m_initial_time << " sec";
        log << t.str();
    }
    log << endreq;

    if( m_tickCount!=0){
        double interval = m_lasttime - m_history->start;
        if( interval<=0) return sc;
        m_history->finish( m_lasttime, interval);
        m_rootTupleSvc->storeRowFlag(this->m_root_tree.value(), true);
    }
    m_history->begin(m_lasttime);

    m_tickCount++;
    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode ExposureAlg::finalize(){
    // finish up
    if( m_tickCount>0 ) execute();
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "Processed " << m_tickCount << " ticks" << endreq;

    delete m_history;
    return sc;
}


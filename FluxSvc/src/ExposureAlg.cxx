/*** @file ExposureAlg.cxx
    @brief declaration and implementation of class ExposureAlg

    $Header$

    note: reverted to 1.44 by THB on 11/10
*/
// Include files
#include "FluxSvc/PointingInfo.h"

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
#include "facilities/Observer.h"

// access to accumulated livetime
#include "Trigger/ILivetimeSvc.h"

#include <cassert>
#include <vector>
#include <fstream>
#include <iomanip>
#include <map>

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


private: 

    double m_lasttime; //time value to hold time between events;
    double m_last_livetime;
    double m_initial_time; 

    StringProperty m_root_tree;
    StringProperty m_pointing_history_input_file;
    StringProperty m_clockName;

    IntegerProperty m_print_frequency;
    IFluxSvc*   m_fluxSvc;

    int         m_tickCount; // number of ticks processed
    INTupleWriterSvc* m_rootTupleSvc;;
    ILivetimeSvc *    m_livetimeSvc;

    PointingInfo m_history;

    // this to support detection of SAA boundary
    bool m_insideSAA;
    ObserverAdapter< ExposureAlg > m_observer; //obsever tag
    int askGPS(); ///< notified when position changes
    
    void createEntry(); // called for tick, or SAA transition

    astro::GPS * m_gps;
    std::map<std::string, int> m_ticks; // ticks by clock name
};
//------------------------------------------------------------------------

static const AlgFactory<ExposureAlg>  Factory;
const IAlgFactory& ExposureAlgFactory = Factory;

//------------------------------------------------------------------------
//! ctor
ExposureAlg::ExposureAlg(const std::string& name, ISvcLocator* pSvcLocator)
: Algorithm(name, pSvcLocator)
, m_lasttime(0)
, m_tickCount(0)
, m_insideSAA(false)
{
    // declare properties with setProperties calls
    declareProperty("root_tree",m_root_tree="pointing_history"); //doesn't work???

    declareProperty("PrintFrequency", m_print_frequency=1);
    declareProperty("clockName"    , m_clockName = "clock" );// the name to check for
}

//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode ExposureAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

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

    m_initial_time =m_lasttime = startTime.value();

    // get a pointer to RootTupleSvc 
    if( (sc = service("RootTupleSvc", m_rootTupleSvc, true) ). isFailure() ) {
            log << MSG::ERROR << " failed to get the RootTupleSvc" << endreq;
            return sc;
    }
    // get a pointer to LivetimeSvc 
    if( (service("LivetimeSvc", m_livetimeSvc, true) ). isFailure() ) {
        log << MSG::ERROR << " LivetimeSvc is not available" << endreq;
        return StatusCode::FAILURE;
    }

    m_history.setFT2Tuple(m_rootTupleSvc, m_root_tree.value());

    log << MSG::INFO << "Using the clock \""<< m_clockName<< "\" to generate FT2 entries" << endreq;

    // attach an observer to be notified when orbital position changes
        // set callback to be notified when the position changes
    m_observer.setAdapter( new ActionAdapter<ExposureAlg>
        (this, &ExposureAlg::askGPS) );

    m_fluxSvc->attachGpsObserver(&m_observer);
    m_gps = astro::GPS::instance();

    return sc;
}
//------------------------------------------------------------------------
int ExposureAlg::askGPS()
{
    // callback from GPS when position has been updated
    astro::EarthCoordinate pos = astro::GPS::instance()->earthpos();
    bool inside = pos.insideSAA();
    if( inside != m_insideSAA){
        // changed: 

        m_insideSAA= inside;
        createEntry();
    }

    return 0; // can't be void in observer pattern
}

//------------------------------------------------------------------------
//! process a tick
StatusCode ExposureAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    IFlux* flux=m_fluxSvc->currentFlux();
        
    std::string name( flux->name());

    log << MSG::DEBUG << "Clock tick at "<< m_gps->time()<<", clock: "<<name << endreq;

    m_ticks[name]++;

    if(m_clockName.value() == name){
        createEntry();
        setFilterPassed(false); // since this is on a branch, and we want the sequence to fail
    } 

    return sc;
}
//------------------------------------------------------------------------
void ExposureAlg::createEntry()
{
    m_lasttime = astro::GPS::instance()->time();
    if( m_tickCount!=0){
        double interval = m_lasttime - m_history.start_time();
        if( interval<=1e-10) return;
        m_history.finish( m_lasttime, m_livetimeSvc->livetime()-m_last_livetime);

        m_rootTupleSvc->saveRow(this->m_root_tree.value());
    }
    m_last_livetime = m_livetimeSvc->livetime();

    // start new entry
    m_history.set();

    if(  m_tickCount% m_print_frequency==0){
        MsgStream   log( msgSvc(), name() );
        log << MSG::INFO;
        if( log.isActive() ){
            std::stringstream t;
            t   << "tick at " << std::setprecision(10)
                << m_lasttime - m_initial_time << " sec"
                << std::setprecision(3)
                << ", lat, lon, magLat, zenithTheta, SAA = " 
                << m_history.lat_geo << ", " 
                << m_history.lon_geo << ", "
                << m_history.lat_mag << ", " 
                << ( fabs(m_history.zenith_scz)<1e-8? 0: m_history.zenith_scz) << ", "
                << m_history.in_saa;
            log << t.str();
        }
        log << endreq;
    }

    m_tickCount++;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode ExposureAlg::finalize(){
    // finish up
    if( m_tickCount>0 ) createEntry();

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    log << MSG::INFO <<" clock Name     ticks"<<endreq;;
    for(std::map<std::string,int>::const_iterator im=m_ticks.begin(); im !=m_ticks.end(); ++im) {
        log<< MSG::INFO 
            <<std::setw(15) << std::left  << im->first
            << std::setw(10)<< std::right << im->second << endreq;
    }
    log << std::endl;

    return sc;
}


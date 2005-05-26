
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

#include <cassert>
#include <vector>
#include <fstream>
#include <iomanip>

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
    double m_initial_time; 

    StringProperty m_root_tree;
    StringProperty m_pointing_history_input_file;

    IFluxSvc*   m_fluxSvc;

    int         m_tickCount; // number of ticks processed
    INTupleWriterSvc* m_rootTupleSvc;;

    PointingInfo m_history;
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

    m_history.setFT2Tuple(m_rootTupleSvc, m_root_tree.value());

    return sc;
}

//------------------------------------------------------------------------
//! process an event
StatusCode ExposureAlg::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    IFlux* flux=m_fluxSvc->currentFlux();
    m_lasttime = flux->time();
        
    std::string particleName = flux->particleName();


    if(particleName != "TimeTick" && particleName != "Clock"){
        return sc;
    }



    if( m_tickCount!=0){
        double interval = m_lasttime - m_history.start_time();
        if( interval<=0) return sc;
        m_history.finish( m_lasttime, interval);
        m_rootTupleSvc->storeRowFlag(this->m_root_tree.value(), true);
    }
    m_history.set(m_lasttime);

    log << MSG::INFO;
    if( log.isActive() ){
        std::stringstream t;
        t   << "tick at " << std::setprecision(10)
            << m_lasttime - m_initial_time << " sec"
            << std::setprecision(3)
            << ", lat, lon, B, L = " 
            << m_history.lat_geo << ", " 
            << m_history.lon_geo << ", "
            << m_history.B << ", " 
            << m_history.L;
        log << t.str();
    }
    log << endreq;

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

    return sc;
}


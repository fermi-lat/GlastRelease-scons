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
#include "Event/MonteCarlo/D2Entry.h"

//flux
#include "FluxSvc.h"
#include "FluxSvc/IFlux.h"
#include "GPS.h"

#include <cassert>
#include <vector>
#include <fstream>

/** 
* \class ExposureAlg
*
* \brief This is an Algorithm designed to get information about LAT
* position, exposure and livetime from FluxSvc and use it to put information onto the TDS about
* LAT pointing and location characteristics, effectively generating the D2 database.  The "TimeCandle" 
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
    StringProperty m_source_name;
    StringProperty m_file_name;

    IFluxSvc*   m_fluxSvc;
    IFlux *     m_flux;

    std::ostream* m_out;  //for output that looks like the stuff from the astro orbit model test.
    int         m_tickCount; // number of ticks processed
        
    
};
//------------------------------------------------------------------------

static const AlgFactory<ExposureAlg>  Factory;
const IAlgFactory& ExposureAlgFactory = Factory;

//------------------------------------------------------------------------
//! ctor
ExposureAlg::ExposureAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator), m_out(0), m_tickCount(0)
{
    // declare properties with setProperties calls
    declareProperty("source_name",  m_source_name="default");
    declareProperty("file_name",  m_file_name="");
    
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

    if(! m_file_name.value().empty() ){
        m_out = new std::ofstream(m_file_name.value().c_str());
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
    
    Event::McParticleCol* pcol = new Event::McParticleCol;
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
        currentTime = m_fluxSvc->currentFlux()->time();
    }   

    //now, only do the rest of this algorithm if we have a timetick particle.
    std::string particleName = m_fluxSvc->currentFlux()->particleName();
    if(particleName != "TimeTick"){
        log << MSG::DEBUG << particleName << " Not a timetick particle, no D2Database entry created, continuing..." << endreq;
        return StatusCode::SUCCESS;
    }

    //by now, we should know that we have the appropriate particle to make a D2Entry with.
    double secondsperday = 60.*60.*24.;

    // here we get the time characteristics

    
    EarthOrbit orb; //for the following line - this should have a better implementation.
    double julianDate = orb.dateFromSeconds(currentTime);

    //NOTE: this gets an interval from the last time that a TimeTick particle came to this one.
    //in other words, the timeTick particles define the beginning and ends of intervals.
    double intrvalstart = m_lasttime;
    //then get the end of the interval (this can be made into whatever).
    double intrvalend = currentTime;
    //get the livetime from the source itself.
    double livetime = intrvalend - m_lasttime;
    //..and reset the time of this event to be the "last time" for next time.
    m_lasttime = intrvalend;

    intrvalstart = orb.dateFromSeconds(intrvalstart);
    intrvalend = orb.dateFromSeconds(intrvalend);

    //and here the pointing characteristics of the LAT.
    GPS::instance()->getPointingCharacteristics(currentTime);
    //EarthOrbit orbt;
    Hep3Vector location = GPS::instance()->position(currentTime);
    
    // hold onto the cartesian location of the LAT
    double posx = location.x(); 
    double posy = location.y(); 
    double posz = location.z(); 
    std::cout << std::endl;
    double rax = GPS::instance()->RAX();
    double raz = GPS::instance()->RAZ();
    double decx = GPS::instance()->DECX();
    double decz = GPS::instance()->DECZ();
    double razenith = GPS::instance()->RAZenith();
    double deczenith = GPS::instance()->DECZenith();
    
    EarthCoordinate earthpos(location,julianDate);
    double lat = earthpos.latitude();
    double lon = earthpos.longitude();
    double alt = earthpos.altitude();
    bool SAA = earthpos.insideSAA();
    
    SolarSystem sstm;
    
    double ramoon = sstm.direction(astro::SolarSystem::Moon,julianDate).ra();
    double decmoon = sstm.direction(astro::SolarSystem::Moon,julianDate).dec();
    double rasun = sstm.direction(astro::SolarSystem::Sun,julianDate).ra();
    double decsun = sstm.direction(astro::SolarSystem::Sun,julianDate).dec();
    
    // Here the TDS is prepared to receive hits vectors
    // Check for the MC branch - it will be created if it is not available
    
    DataObject *mc = new Event::D2EntryCol;
    sc=eventSvc()->registerObject(EventModel::MC::Event , mc);
    // THB: why is this commented out?
    //if(sc.isFailure()) {
    //    log << MSG::ERROR << EventModel::MC::Event  <<" could not be registered on data store" << endreq;
    //    return sc;
    //}
    
    // Here the TDS receives the exposure data
    Event::D2EntryCol* exposureDBase = new Event::D2EntryCol;
    sc=eventSvc()->registerObject(EventModel::MC::D2EntryCol , exposureDBase);
    if(sc.isFailure()) {
        log << MSG::ERROR << EventModel::MC::D2EntryCol  <<" could not be entered into existing data store" << endreq;
        return sc;
    }
    
    Event::D2Entry* entry = new Event::D2Entry;
    
    exposureDBase->push_back(entry);

    entry->init(posx, posy, posz,rax,raz,decx,decz,razenith,deczenith,lat,
        lon,alt,intrvalstart,intrvalend,
        livetime,ramoon,decmoon,rasun,decsun,
        SAA);
    
    // now we'll retreive the data from the TDS as a check.
    Event::D2EntryCol* elist = new Event::D2EntryCol;
    eventSvc()->retrieveObject("/Event/MC/D2EntryCol",(DataObject *&)elist);
    
    Event::D2EntryCol::iterator curEntry = (*elist).begin();
    //some test output - to show that the data got onto the TDS
    (*curEntry)->writeOut(log);

    SkyDir curDir(raz,decz);
    SkyDir xDir(rax,decx);

    SkyDir sunDir(rasun,decsun);
    //Rotation galtoglast(m_fluxSvc->transformGlastToGalactic(currentTime).inverse);
    sunDir()=(m_fluxSvc->transformGlastToGalactic(currentTime).inverse())*sunDir();
     
    //and here's the file output.
    if( m_out !=0) {
        std::ostream& out = *m_out;
        out<<intrvalstart <<'\t';
        out<<intrvalend <<'\t';
        out<< posx<<'\t';
        out<< posy<<'\t';
        out<< posz<<'\t';
        out<<raz<<'\t';
        out<<decz<<'\t';
        out<<rax<<'\t';
        out<<decx<<'\t';
        out<<razenith<<"\t";
        out<<deczenith <<'\t';
        out<<"1"<<'\t';
        out<<livetime <<'\t';
        out<<SAA<<'\t';
        out<<lon <<'\t';
        out<<lat <<'\t';
        out<<alt <<'\t';
        out<<curDir.l() <<'\t';
        out<<curDir.b()<<'\t';
        out<<xDir.l() <<'\t';
        out<<xDir.b()<<'\t';
        out<< rasun <<"\t"<<decsun  <<'\t';
        out<< sunDir().x() <<"\t"<< sunDir().y()  <<"\t"<< sunDir().z()  <<'\t';
        out<<ramoon <<"\t"<<decmoon   <<std::endl;
        
        
    }
        setFilterPassed( false );
    log << MSG::DEBUG << "ExposureAlg found a TimeTick particle, ended this execution after making a record, filterpassed = " << filterPassed() << endreq;

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


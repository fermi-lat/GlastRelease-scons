// $Header$

// Include files
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"
#include "GaudiKernel/SmartRefVector.h"

#include "astro/SkyDir.h"
#include "astro/EarthCoordinate.h"
#include "astro/SolarSystem.h"
// Event for creating the McEvent stuff
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/D2Entry.h"

//flux
#include "FluxSvc.h"
#include "FluxSvc/IFlux.h"
#include "Spectrum.h"
#include "SpectrumFactory.h"
#include "GPS.h"

//for instantiation of the new spectrum.
#include "FluxSvc/ISpectrum.h"
#include "FluxSvc/ISpectrumFactory.h" 

#include "astro/EarthOrbit.h"

#include "CLHEP/Vector/LorentzVector.h"

#include <cassert>
#include <vector>

#include "ExposureAlg.h"
//------------------------------------------------------------------------


static const AlgFactory<ExposureAlg>  Factory;
const IAlgFactory& ExposureAlgFactory = Factory;

static std::ofstream m_out("orbitFromAlg.out");
//------------------------------------------------------------------------
//! ctor
ExposureAlg::ExposureAlg(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
{
    // declare properties with setProperties calls
    declareProperty("source_name",  m_source_name="default");
    
}
//-------------------------------------------------------------------------
//! simple prototype of a spectrum class definition - this implememts all of the necessary methods in an ISpectrum
/*class TimeCandle : public ISpectrum {
public:
    TimeCandle(const std::string& params){};
    
    virtual const char * particleName()const{return "gamma";}
    
    /// return a title describing the spectrum	
    virtual std::string title()const{ return "standard time candle";}
    
    /// calculate flux 
    //virtual double flux() const { return 1.0;}
    
    virtual double flux(double time) const {return 1.0;}
    
    /// return kinetic energy (in GeV) 
    /// @param r Number between 0 and 1 for min to max energy
    virtual float operator()(float r)const{
        return r;
    }
    
    virtual std::pair<float,float> dir(float)const{
        return std::make_pair<float,float>(1.0,0.0);
    }     
    
    virtual std::pair<double,double> dir(double energy, HepRandomEngine* engine){return std::make_pair<double,double>(1.0,0.0);}
    
    double energySrc(HepRandomEngine* engine, double time){  return (*this)(engine->flat());}
    
    double interval (double time)
    {        
        return 5.;
    }
    
};*/
//! set this spectrum up to be used.
void ExposureAlg::makeTimeCandle(IFluxSvc* fsvc){
    //static RemoteSpectrumFactory<TimeCandle> timecandle(fsvc);
}

//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode ExposureAlg::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;

    //m_out = new std::ofstream("orbitFromAlg.out");
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    if ( service("FluxSvc", m_fluxSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the FluxSvc!" << endreq;
        return StatusCode::FAILURE;
    }
    
    // set up the standard time candle spectrum
    makeTimeCandle(m_fluxSvc);
    
    log << MSG::INFO << "loading source..." << endreq;
    
    sc =  m_fluxSvc->source(m_source_name, m_flux);
    if( sc.isFailure()) {
        log << MSG::ERROR << "Could not find flux " << m_source_name << endreq;
        return sc;
    }
    log << MSG::INFO << "Source: "<< m_flux->title() << endreq;
    
    log << MSG::INFO << "Source title: " << m_flux->title() << endreq;
    log << MSG::INFO << "        area: " << m_flux->targetArea() << endreq;
    log << MSG::INFO << "        rate: " << m_flux->rate() << endreq;
    
    if ( service("ParticlePropertySvc", m_partSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the ParticlePropertySvc!" << endreq;
        return StatusCode::FAILURE;
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

    //NOTE: this gets an interval from the last time that a TimeTick particle came to this one.
    //in other words, the timeTick particles define the beginning and ends of intervals.
    double intrvalstart = m_lasttime;
    //then get the end of the interval (this can be made into whatever).
    double intrvalend = currentTime;
    //get the livetime from the source itself.
    double livetime = intrvalend - m_lasttime;
    //..and reset the time of this event to be the "last time" for next time.
    m_lasttime = intrvalend;
    
    //and here the pointing characteristics of the LAT.
    GPS::instance()->getPointingCharacteristics(/*location,20*/currentTime/secondsperday);
    //EarthOrbit orbt;
    Hep3Vector location = GPS::instance()->position();//( orbt.position(currentTime/secondsperday) );
    
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
    
    EarthCoordinate earthpos(location,currentTime/secondsperday);
    double lat = earthpos.latitude();
    double lon = earthpos.longitude();
    double alt = earthpos.altitude();
    bool SAA = earthpos.insideSAA();
    
    SolarSystem sstm;
    
    double ramoon = sstm.direction(astro::SolarSystem::Moon,currentTime/secondsperday).ra();
    double decmoon = sstm.direction(astro::SolarSystem::Moon,currentTime/secondsperday).dec();
    double rasun = sstm.direction(astro::SolarSystem::Sun,currentTime/secondsperday).ra();
    double decsun = sstm.direction(astro::SolarSystem::Sun,currentTime/secondsperday).dec();
    
    // Here the TDS is prepared to receive hits vectors
    // Check for the MC branch - it will be created if it is not available
    
    DataObject *mc = new Event::D2EntryCol;
    //eventSvc()->retrieveObject("/Event/MC", mc);
    sc=eventSvc()->registerObject(EventModel::MC::Event , mc);
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
    
    entry->init(posx, posy, posz,rax,raz,decx,decz,razenith,deczenith,lat,
        lon,alt,intrvalstart,intrvalend,
        livetime,ramoon,decmoon,rasun,decsun,
        SAA);
    
    exposureDBase->push_back(entry);
    
    // now we'll retreive the data from the TDS as a check.
    Event::D2EntryCol* elist = new Event::D2EntryCol;
    eventSvc()->retrieveObject("/Event/MC/D2EntryCol",(DataObject *&)elist);
    
    Event::D2EntryCol::iterator curEntry = (*elist).begin();
    //some test output - to show that the data got onto the TDS
    (*curEntry)->writeOut(log);
//WRITEOUT------------------------------------------------------------------------------------------

     
        //m_out.flags(ios::fixed);
        m_out<<intrvalstart <<'\t';
        m_out<<intrvalend <<'\t';
        m_out<< posx<<'\t';
        m_out<< posy<<'\t';
        m_out<< posz<<'\t';
        m_out<<raz<<'\t';
        m_out<<decz<<'\t';
        m_out<<rax<<'\t';
        m_out<<decx<<'\t';
        m_out<<razenith<<"\t";
        m_out<<deczenith <<'\t';
        m_out<<"1"<<'\t';
        m_out<<livetime <<'\t';
        m_out<<SAA<<'\t';
        m_out<<lon <<'\t';
        m_out<<lat <<'\t';
        m_out<<alt <<'\t';
        m_out<< rasun <<"\t"<<decsun  <<'\t';
        m_out<<ramoon <<"\t"<<decmoon   <<std::endl;
                                
                                

//WRITEOUT------------------------------------------------------------------------------------------------

    setFilterPassed( false );
    log << MSG::DEBUG << "ExposureAlg found a TimeTick particle, ended this execution after making a record, filterpassed = " << filterPassed() << endreq;

    
    return sc/*return StatusCode::STOP*/;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode ExposureAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    m_out.close();
    
    return sc;
}


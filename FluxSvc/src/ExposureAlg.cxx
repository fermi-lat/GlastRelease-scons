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
class TimeCandle : public ISpectrum {
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
    
};
//! set this spectrum up to be used.
void ExposureAlg::makeTimeCandle(IFluxSvc* fsvc){
    static RemoteSpectrumFactory<TimeCandle> timecandle(fsvc);
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
    
    if(m_fluxSvc->currentFlux() == m_flux){
        m_flux->generate();
    }else{
        m_flux = m_fluxSvc->currentFlux();
        m_flux->generate();
    }
    
    double secondsperday = 60.*60.*24.;
    EarthOrbit orbt;
    Hep3Vector location( orbt.position(m_flux->time()/secondsperday) );
    
    // hold onto the cartesian location of the LAT
    double posx = location.x(); 
    double posy = location.y(); 
    double posz = location.z(); 
    
    // here we get the time characteristics
    double intrvalstart = m_lasttime;
    //then get the end of the interval (this can be made into whatever).
    double intrvalend = m_flux->time();
    //get the livetime from the source itself.
    double livetime = intrvalend - m_lasttime;
    //..and reset the time of this event to be the "last time" for next time.
    m_lasttime = intrvalend;
    
    //and here the pointing characteristics of the LAT.
    GPS::instance()->getPointingCharacteristics(location,20);
    std::cout << std::endl;
    double rax = GPS::instance()->RAX();
    double raz = GPS::instance()->RAZ();
    double decx = GPS::instance()->DECX();
    double decz = GPS::instance()->DECZ();
    double razenith = GPS::instance()->RAZenith();
    double deczenith = GPS::instance()->DECZenith();
    
    EarthCoordinate earthpos(location,m_flux->time());
    double lat = earthpos.latitude();
    double lon = earthpos.longitude();
    double alt = earthpos.altitude();
    bool SAA = earthpos.insideSAA();
    
    SolarSystem sstm;
    
    double ramoon = sstm.direction(astro::SolarSystem::Moon,m_flux->time()).ra();
    double decmoon = sstm.direction(astro::SolarSystem::Moon,m_flux->time()).dec();
    double rasun = sstm.direction(astro::SolarSystem::Sun,m_flux->time()).ra();
    double decsun = sstm.direction(astro::SolarSystem::Sun,m_flux->time()).dec();
    
    // Here the TDS is prepared to receive hits vectors
    // Check for the MC branch - it will be created if it is not available
    
    DataObject *mc = new Event::D2EntryCol;
    //eventSvc()->retrieveObject("/Event/MC", mc);
    sc=eventSvc()->registerObject(EventModel::MC::Event , mc);
    if(sc.isFailure()) {
        log << MSG::ERROR << EventModel::MC::Event  <<" could not be registered on data store" << endreq;
        return sc;
    }
    
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
    
    return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode ExposureAlg::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    return sc;
}


// $Header$

// Include files
#include "FluxSvc/IFluxSvc.h"
#include "FluxSvc/IFlux.h"


// Event for creating the McEvent stuff
//#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"


// Gaudi system includes
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include <list>
#include <string>
#include <vector>
#include "GaudiKernel/ParticleProperty.h"

//#include "FluxAlg.h"
/*! \class FluxTestAlg
\brief 
In addition to the normal Gaudi JobOptions requirements, there are:
FluxSvc.source_lib, which should contain the relevant xml files to be used. and
FluxTestAlg.source_name, which holds the name of the desired spectrum.
*/
//class ParticleProperty;

class FluxTestAlg : public Algorithm {
    
public:
    //! Constructor of this form must be provided
    FluxTestAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
    
private:
    IFlux* m_flux;
    std::string m_source_name;
    IParticlePropertySvc * m_partSvc;
    std::ostream* m_out;  //for output that looks like the stuff from the astro orbit model test.
    std::ostream* m_diffsources; //for output concerning the source used (extra output for diffuse model)

};


static const AlgFactory<FluxTestAlg>  Factory;
const IAlgFactory& FluxTestAlgFactory = Factory;
/*
void WARNING (const char * text ){  std::cerr << "WARNING: " << text << '\n';}
void FATAL(const char* s){std::cerr << "\nERROR: "<< s;}
*/

//------------------------------------------------------------------------------
//
FluxTestAlg::FluxTestAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator){
    
    declareProperty("source_name", m_source_name="default");
}


//------------------------------------------------------------------------------
/*! */
StatusCode FluxTestAlg::initialize() {
    
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initializing..." << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    if ( service("ParticlePropertySvc", m_partSvc).isFailure() ){
            log << MSG::ERROR << "Couldn't find the ParticlePropertySvc!" << endreq;
            return StatusCode::FAILURE;
    }

#if 0 // this does not make sense for testing FluxAlg
    //set the output file.
    m_out = new std::ofstream("TestOutputData.out");
    m_diffsources = new std::ofstream("SourceCharacteristics.out");

    // get the pointer to the flux Service 
    IFluxSvc* fsvc;
    
    // get the service
    StatusCode sc = service("FluxSvc", fsvc);
    
    if( sc.isFailure()) {
        log << MSG::ERROR << "Could not find FluxSvc" << endreq;
        return sc;
    }
    
    //uncomment this to set the rocking method.
    //fsvc->setRockType(GPS::UPDOWN);
    
    log << MSG::INFO << "loading source..." << endreq;
    
    
    sc =  fsvc->source(m_source_name, m_flux);
    if( sc.isFailure()) {
        log << MSG::ERROR << "Could not find flux " << m_source_name << endreq;
        return sc;
    }
    
    // then do the output here.
    log << MSG::INFO << "start of other loops" << endreq;
    log << MSG::INFO << "Source title: " << m_flux->title() << endreq;
    log << MSG::INFO << "       area: " << m_flux->targetArea() << endreq;
    log << MSG::INFO << "       rate: " << m_flux->rate() << endreq;
    
#endif   
    
    return StatusCode::SUCCESS;
}


//------------------------------------------------------------------------------
StatusCode FluxTestAlg::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );    
    
    
    Event::McParticleCol* pcol = 0;
    eventSvc()->retrieveObject(EventModel::MC::McParticleCol, (DataObject *&)pcol);
    
    if(pcol==0){ 
        log << MSG::ERROR << " Could not find "<< EventModel::MC::McParticleCol << endreq;
            return StatusCode::FAILURE;
    }

    SmartDataPtr<Event::EventHeader> header(eventSvc(), EventModel::EventHeader);
    if( 0==header) {  
        log << MSG::ERROR << " Could not find "<< EventModel::EventHeader << endreq;
            return StatusCode::FAILURE;
    }



    Event::McParticleCol::iterator elem = (*pcol).begin();
    HepVector3D d = (*elem)->initialFourMomentum().v().unit();
    HepVector3D p = (*elem)->finalPosition();
    
    double energy = (*elem)->initialFourMomentum().e();
    /*StdHepId*/int pID = (*elem)->particleProperty();
    std::string partName = m_partSvc->findByStdHepID(pID)->particle();
    
    log << MSG::INFO << partName
        << "(" << energy
        << " MeV), Launch: " 
        << "(" << p.x() <<", "<< p.y() <<", "<<p.z()<<")" 
        << " Dir " 
        << "(" << d.x() <<", "<< d.y() <<", "<<d.z()<<")"
         << ",  Elapsed Time = " << header->time()
        << endreq;

#if 0 // enable for tests: see above
    double theta = acos(d.z()/(d.mag()));

    double phi = atan2(d.y(),d.x());

    //and here's the file output.
    std::ostream& out = *m_out;
        out<<m_flux->time() <<'\t';
        out<<d.x() <<'\t';
        out<< d.y()<<'\t';
        out<< d.z()<<'\t';
        out<< theta<<'\t';
        out<< phi<<'\t' << std::endl;
#endif
    return sc;
}


//------------------------------------------------------------------------------
StatusCode FluxTestAlg::finalize() {
    delete m_out;
    std::ostream& diffsources = *m_diffsources;
    m_flux->writeSourceCharacteristic(diffsources);
    delete m_diffsources;
    return StatusCode::SUCCESS;
}








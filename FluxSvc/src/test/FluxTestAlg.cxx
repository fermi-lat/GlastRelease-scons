// $Header$

// Include files
#include "FluxSvc/IFluxSvc.h"
#include "FluxSvc/IFlux.h"


// Event for creating the McEvent stuff
//#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/EventModel.h"

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


#include "GalacticHist.h"

//#include "FluxAlg.h"
/*! \class FluxTestAlg
\brief 

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
    
    typedef struct{
        int x;
        int y;
        double amount;
    }exposureSet;
    
    IFlux* m_flux;
    std::string m_source_name;
    IParticlePropertySvc * m_partSvc;
    double m_exposedArea[180][90];
    double m_currentTime;
    double m_passedTime;  //time passed during this event
    std::vector<exposureSet> findExposed(double l,double b, double deltat);
    void addToTotalExposure(std::vector<exposureSet>);
    void displayExposure();
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
    
    // get the pointer to the flux Service 
    IFluxSvc* fsvc;
    
    // get the service
    StatusCode sc = service("FluxSvc", fsvc);
    
    if( sc.isFailure()) {
        log << MSG::ERROR << "Could not find FluxSvc" << endreq;
        return sc;
    }
    
    
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
    
    
    
    return sc;
}


//------------------------------------------------------------------------------
StatusCode FluxTestAlg::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );    
    
    //Event::McParticleCol*  pcol2= SmartDataPtr<Event::McParticleCol>(eventSvc(), "/Event/MC/McParticleCol");
    
    Event::McParticleCol* pcol = new Event::McParticleCol;
    eventSvc()->retrieveObject("/Event/MC/McParticleCol",(DataObject *&)pcol);
    
    HepVector3D p,d;
    double energy;
    std::string partName;
    //only make a new source if one does not already exist.
    if(pcol==0){
        m_flux->generate();
        p = m_flux->launchPoint()*10.;
        d = m_flux->launchDir();
        energy = m_flux->energy();
        partName = m_flux->particleName();
    }else{
        Event::McParticleCol::iterator elem = (*pcol).begin();
        d = (*elem)->initialFourMomentum().v()*10;
        p = (*elem)->finalPosition();
        
        energy = (*elem)->initialFourMomentum().e();
        /*StdHepId*/int pID = (*elem)->particleProperty();
        if ( service("ParticlePropertySvc", m_partSvc).isFailure() ){
            log << MSG::ERROR << "Couldn't find the ParticlePropertySvc!" << endreq;
            return StatusCode::FAILURE;
        }
        partName = m_partSvc->findByStdHepID(pID)->particle();
    }
    
    /*
    log << MSG::INFO << partName
    << "(" << energy
    << " GeV), Launch: " 
    << "(" << p.x() <<", "<< p.y() <<", "<<p.z()<<")" 
    << " Dir " 
    << "(" << d.x() <<", "<< d.y() <<", "<<d.z()<<")"
    // << ",  Elapsed Time = " << m_flux->time()
    << endreq;   */
    
    // get the pointer to the flux Service 
    IFluxSvc* fsvc;
    // get the service
    sc = service("FluxSvc", fsvc);
    if( sc.isFailure()) {
        log << MSG::ERROR << "Could not find FluxSvc" << endreq;
        return sc;
    }
    
    HepVector3D pointingin = d;//(0,0,1);
    pointingin = (fsvc->transformGlastToGalactic(m_flux->time()))*pointingin;
    
    //log << MSG::INFO
    //        << "(" << pointingin.x() <<", "<< pointingin.y() <<", "<<pointingin.z()<<")" 
    //        << endreq;
    
    //we want to make this into l and b now.
    double l,b;
    l = atan(pointingin.x()/pointingin.z());
    b = atan(pointingin.y()/pointingin.z());
    
    l *= 360./M_2PI;
    b *= 360./M_2PI;
    
    l+= 180;
    b+= 90;
    
    log << MSG::INFO
        << "(" << "l = " << l << ", b = " << b <<")" 
        << endreq;

    m_passedTime = (m_flux->time())-m_currentTime;
    m_currentTime = m_flux->time();

    std::vector<exposureSet> exposed;
    exposed = findExposed(l,b,m_passedTime);
    addToTotalExposure(exposed);
    
    
    //m_flux->pass(10.);
    return sc;
}


//------------------------------------------------------------------------------
StatusCode FluxTestAlg::finalize() {
    displayExposure();
    return StatusCode::SUCCESS;
}



std::vector<FluxTestAlg::exposureSet> FluxTestAlg::findExposed(double l,double b,double deltat){
    std::vector<exposureSet> returned;
    double angularRadius = 10;
    for(int i= l-angularRadius ; i<=l+angularRadius ; i+=2){
        for(int j= b-angularRadius ; j<=b+angularRadius ; j+=2){
            
            if((pow(l-i,2)+pow(b-j,2) <= pow(angularRadius,2))){
                //set up the point, and stick it into the vector
                exposureSet point;
                point.x = i/2; //yes, this is doing an implicit cast.
                point.y = j/2;  //these should be divided by two, but they're being shrunk for the current display.
                point.amount = deltat;
                returned.push_back(point);
            }
        }
    }
    return returned;
}

void FluxTestAlg::addToTotalExposure(std::vector<FluxTestAlg::exposureSet> toBeAdded){
    std::vector<exposureSet>::iterator iter = toBeAdded.begin();
    int x,y;
    double amount;
    if(toBeAdded.size()){
        for( ; iter!=toBeAdded.end() ; iter++){
            x = (*iter).x;
            y = (*iter).y;
            amount = (*iter).amount;
            m_exposedArea[x][y] += amount;
            
        }
    }else{
        
        std::cout << "error in addToTotalExposure - null vector input" <<std::endl;
    }
}

void FluxTestAlg::displayExposure(){
    //make the file
    std::ofstream out_file("data.dat", std::ios::ate);

    int i,j;
    for(i=0 ; i<180 ; i++){
        std::strstream out;
        for(j=0 ; j<90 ; j++){
            out << m_exposedArea[i][j] << " ";
            out_file << i << "  " << j << "  " << m_exposedArea[i][j] << std::endl;
        }
        out << std::endl;
        std::cout << out.str();
    }
    out_file.close();
    //then use galacticHist to display the stuf in the file
    galacticHist abc;
    abc.test();
}




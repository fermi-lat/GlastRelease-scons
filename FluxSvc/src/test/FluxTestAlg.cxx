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
#include <math.h>
#include "GaudiKernel/ParticleProperty.h"

#include "FluxSvc/ISpectrum.h"
#include "FluxSvc/ISpectrumFactory.h" 

#include "CLHEP/Random/RandFlat.h"


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
    double energy;
    double m_glastExposureAngle;
    int m_exposureFunction;
    double m_pointSpread;
    int m_projectionType;
    int m_rockingMode;
    
    double m_exposedArea[360][180];
    double m_angleExposure[100];
    double m_currentTime;
    double m_passedTime;  //time passed during this event
    void findExposed(double l,double b, double deltat, Rotation glastToGal);
    void collectZenithHist(double l,double b, double deltat);
    void displayExposure();
    void rootDisplay();
    void zenithHistDisplay();
    void makeTimeCandle(IFluxSvc* fsvc);
    double zenithExposure(double zenithAngle);
    double energyExposure(double energy);
    double azimuthExposure(double azimuth);
    
    /// transform linear l,b coordinates into transformed equal-area coords.
    std::pair<double,double> hammerAitoff(double l,double b);
    
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
    declareProperty("glastExposureAngle", m_glastExposureAngle=10.);
    declareProperty("exposureFunction", m_exposureFunction=0);
    declareProperty("pointSpread", m_pointSpread=0.);
    declareProperty("projectionType", m_projectionType=0);
    declareProperty("rockingMode", m_rockingMode=0);
    
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
    
    // set up the standard time candle spectrum
    makeTimeCandle(fsvc);    
    
    log << MSG::INFO << "loading source..." << endreq;
    
    
    sc =  fsvc->source(m_source_name, m_flux);
    if( sc.isFailure()) {
        log << MSG::ERROR << "Could not find flux " << m_source_name << endreq;
        return sc;
    }

    //for testing - offset of start of data run
    //m_flux->pass(30.);

    // then do the output here.
    log << MSG::INFO << "start of other loops" << endreq;
    log << MSG::INFO << "Source title: " << m_flux->title() << endreq;
    log << MSG::INFO << "       area: " << m_flux->targetArea() << endreq;
    log << MSG::INFO << "       rate: " << m_flux->rate() << endreq;
    
    return sc;
}


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
/// set this spectrum up to be used.
void FluxTestAlg::makeTimeCandle(IFluxSvc* fsvc){
    static RemoteSpectrumFactory<TimeCandle> timecandle(fsvc);
}

//------------------------------------------------------------------------------
StatusCode FluxTestAlg::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );    
    
    //Event::McParticleCol*  pcol2= SmartDataPtr<Event::McParticleCol>(eventSvc(), "/Event/MC/McParticleCol");
    
    Event::McParticleCol* pcol = new Event::McParticleCol;
    eventSvc()->retrieveObject("/Event/MC/McParticleCol",(DataObject *&)pcol);
    
    // get the pointer to the flux Service 
    IFluxSvc* fsvc;
    // get the service
    sc = service("FluxSvc", fsvc);
    if( sc.isFailure()) {
        log << MSG::ERROR << "Could not find FluxSvc" << endreq;
        return sc;
    }
    
    HepVector3D p,d;
    //double energy;
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
    
    
    /* log << MSG::INFO << partName
    << "(" << energy
    << " GeV), Launch: " 
    << "(" << p.x() <<", "<< p.y() <<", "<<p.z()<<")" 
    << " Dir " 
    << "(" << d.x() <<", "<< d.y() <<", "<<d.z()<<")"
    // << ",  Elapsed Time = " << m_flux->time()
    << endreq;   
    */
    
    HepVector3D pointingin(0,0,-1);
    if(m_rockingMode==1){
        /// FOR ROCKING ONE DIRECTION ON EACH ORBIT
        //this needs to be fixed - no hardwired periods!.
        if((int)m_flux->gpsTime()%180 <= 97){
            // Set Rocking angles (zenith rotation, off-zenith)
            fsvc->setOrientation(std::make_pair<double,double>(0.0,30.*M_PI/180.));
        }else{
            fsvc->setOrientation(std::make_pair<double,double>(0.0,-30.*M_PI/180.));
        }
    }else if(m_rockingMode==2){
        /// FOR ROCKING BOTH DIRECTIONS ON EACH ORBIT
        //this needs to be fixed - no hardwired periods!.
        if((int)m_flux->gpsTime()%90 <= 48){
            // Set Rocking angles (zenith rotation, off-zenith)
            fsvc->setOrientation(std::make_pair<double,double>(0.0,30.*M_PI/180.));
        }else{
            fsvc->setOrientation(std::make_pair<double,double>(0.0,-30.*M_PI/180.));
        }
    }else{
        //if we're not using rocking, it makes sense we may be sky mapping:
        pointingin = d;
    }
    
    Rotation glastToGal = fsvc->transformGlastToGalactic(m_flux->gpsTime());
    
    pointingin = (fsvc->transformGlastToGalactic(m_flux->gpsTime()))*pointingin;
    
    //log << MSG::INFO
    //        << "(" << pointingin.x() <<", "<< pointingin.y() <<", "<<pointingin.z()<<")" 
    //        << endreq;
    
    //we want to make this into l and b now.
    double l,b;
    l = atan(pointingin.x()/pointingin.z());
    b = atan(pointingin.y()/pointingin.z());
    //std::cout << "pointingin z is " << pointingin.z() <<std::endl;
    
    //l = asin(pointingin.x());
    b = asin(pointingin.y());
    
    l *= 360./M_2PI;
    b *= 360./M_2PI;
    
    //a serious kluge - this part needs further examination
    if(pointingin.z()<0){
        if(pointingin.x()>=0){
            l=180.+l;
        }else if(pointingin.x()<=0){
            l=-180.+l;
        }
        // l *= -1.;
    }
    
    l+= 180;
    b+= 90;
    
    //log << MSG::INFO
    //    << "(" << "l = " << l << ", b = " << b <<")" 
    //    << endreq;
    
    m_passedTime = (m_flux->time())-m_currentTime;
    m_currentTime = m_flux->time();
    
    if(1){
    findExposed(l,b,m_passedTime,glastToGal);
    }else{
        collectZenithHist(l,b,m_passedTime);
    }
    return sc;
}


//------------------------------------------------------------------------------
StatusCode FluxTestAlg::finalize() {
    if(1){displayExposure();
    }else{
    zenithHistDisplay();
    }
    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
void FluxTestAlg::collectZenithHist(double l,double b, double deltat){
    double watchedl = -15+180;    
    double watchedb = -0+90.;  
    double degToRad = M_PI/180;
    int angularDistance=acos(sin((watchedb-90.)*degToRad)*sin((b-90.)*degToRad)+cos((b-90.)*degToRad)*cos((watchedb-90.)*degToRad)*cos((l-watchedl)*degToRad))/degToRad;
    m_angleExposure[(int)(cos(angularDistance*degToRad)*100)]++;
    return;
}

//------------------------------------------------------------------------------

void FluxTestAlg::findExposed(double l,double b,double deltat, Rotation glastToGal){    
    MsgStream log(msgSvc(), name());
    
    std::vector<exposureSet> returned;
    double angularRadius = m_glastExposureAngle/2.;
    
    //if(m_exposureMode == 6){
    double startl=0; double startb=0;
    double endl=360; double endb=180;
    //now, if you can afford to be more economical when searching here, do so:
    if(m_exposureFunction == 0){
        startl=l;
        startb=b;
        endl=l;
        endb=b;
    }else if(m_exposureFunction == 1 || m_exposureFunction == 2){
        startl=l-angularRadius;
        startb=b-angularRadius;
        endl=l+angularRadius+1;
        endb=b+angularRadius+1;
    }
    //and now start looping over the important area
    for(int i= startl; i<endl ; i++){
        for(int j= startb ; j<endb ; j++){
            
            //if((pow(l-i,2)+pow(b-j,2) <= pow(angularRadius,2))){
            //set up the point, and stick it into the vector
            exposureSet point;// = new exposureSet;
            float correctedl = i;   // Fold into the range [0, 360)
            float correctedb = j;   // Fold into the range [0, 360)                   
            if(correctedl < 0)correctedl+=360;
            if(correctedb < /*-90*/0)correctedb+=180;
            if(correctedl > 360)correctedl-=360;
            if(correctedb > /*90*/180)correctedb-=180;
            
            // determine if the current area is within the exposed area, 
            // and increment it by a number to be multiplied by the integrated time.
            int validpoint=0;
            double exposure;
            if(m_exposureFunction == 0){
                point.x = l; //yes, this is doing an implicit cast.
                point.y = b;
                exposure = 1.;
                validpoint++;
            }else if(m_exposureFunction == 1){
                if((pow(l-i,2)+pow(b-j,2) <= pow(angularRadius,2))){
                    point.x=correctedl;
                    point.y=correctedb;
                    exposure=1.;
                    validpoint++;
                }
            }else if(m_exposureFunction == 2){
                if((pow(l-i,2)+pow(b-j,2) <= pow(angularRadius,2))){
                    
                    point.x = correctedl; //yes, this is doing an implicit cast.
                    point.y = correctedb;
                    
                    //here, the amount depends on the time, and a gaussian in azimuth
                    exposure = pow(2.71828, -(pow(l-i,2)+pow(b-j,2))/(2*pow(angularRadius,2)));
                    //point.amount = gaussian;
                    validpoint++;
                }
            }else if(m_exposureFunction == 3){
                double degToRad = M_PI/180;
                
                if(1.5*acos(sin((j-90.)*degToRad)*sin((b-90.)*degToRad)+cos((b-90.)*degToRad)*cos((j-90.)*degToRad)*cos((l-i)*degToRad))/degToRad  <= pow(angularRadius,1)){
                    point.x=i;
                    point.y=j;
                    exposure=1;
                    validpoint++;
                }
            }else if(m_exposureFunction == 4){
                double degToRad = M_PI/180;
                double azimuth=0;
                double observedl=i-180.;double observedb=j-90.;
                //first, we need to generate a glast local vector corresponding to 
                //the current "observed" l,b point:
                double theta=sqrt(pow(observedb,2)+pow(observedl,2))*M_2PI/360.;
                if (observedl==0.){observedl+=0.000000000001;}  //to fix divide-by-zero errors
                double phi = atan(observedb/observedl);
                //then make the galactic cartesian vector for the "observed" location
                HepVector3D observedgal(-sin(observedl*M_2PI/360.)*cos(observedb*M_2PI/360.) , -sin(observedb*M_2PI/360.) , -cos(observedl*M_2PI/360.)*cos(observedb*M_2PI/360.));
                //and do the transform so we have a "glast-local" direction
                observedgal = (glastToGal.inverse())*observedgal;
                azimuth = atan(observedgal.y()/observedgal.x());
                azimuth*=180./M_PI;
                //if(observedgal.x()>0)azimuth=0;
                    if(observedgal.x()<0 && observedgal.y()<0)azimuth=-180+azimuth;
                    if(observedgal.x()<0 && observedgal.y()>0)azimuth=180+azimuth;
                    //now it is in the range (-180,180)
                //azimuth+=180; //should now be in the range (0,360)
                //place the point to be observed on the correct coordinates                
                point.x=i;
                point.y=j;
                double zenithAngle = acos(sin((j-90.)*degToRad)*sin((b-90.)*degToRad)+cos((b-90.)*degToRad)*cos((j-90.)*degToRad)*cos((l-i)*degToRad))/degToRad;
                //double azimuth=0.;
                exposure=zenithExposure(zenithAngle)*energyExposure(energy)*azimuthExposure(azimuth);
                if(exposure)validpoint++;
            }
            
            
            // Calculate the shift from the point spread function
            double lshift=0,bshift=0,theta,phi;
            if(m_pointSpread){
                theta=m_pointSpread*log10(10.0/(RandFlat::shoot(1.0)));
                phi=RandFlat::shoot(1.0)*M_2PI;
                lshift=theta*sin(phi);
                bshift=theta*cos(phi);
            }
            
            
            /// apply a projection, if desired.
            if(m_projectionType){
                std::pair<double,double> abc = hammerAitoff(correctedl-180.,correctedb-90.);
                point.x = abc.first+lshift; //yes, this is doing an implicit cast.
                point.y = abc.second+bshift;
                point.amount = exposure*deltat*pow(cos((point.y-90.)*M_PI/180.),0.75);
            }else{
                point.x = correctedl+lshift; //yes, this is doing an implicit cast.
                point.y = correctedb+bshift;
                point.amount = exposure*deltat;
            }
            if(validpoint)m_exposedArea[point.x][point.y] += point.amount;
        }
    }
    return;
}

//------------------------------------------------------------
void FluxTestAlg::zenithHistDisplay(){
//make the file
    std::ofstream out_file1("data2.dat", std::ios::ate);
    
    int i;
    for(i=0 ; i<100 ; i++){

            out_file1 << i << "  " << m_angleExposure[i] << std::endl;
    }
    out_file1.close();
    //then use rootDisplay to display the stuff in the file
    std::ofstream out_file("graph2.cxx", std::ios::app);    
    out_file.clear();
    
    out_file << 
        "{\n"
        
        "  FILE *fp;\n"
        "  float ptmp,p[20];\n"
        "  int i, iline=0;\n"; 
    
    out_file <<
        "  TH1D *hist2 = new TH1D(" << '"' << "hist2" << '"' << "," << '"' << "Total Exposure" << '"' << ",100,0.,100.);\n";
    
    out_file <<
        "  hist2->SetXTitle(" << '"' << "cos(theta)(x100)" << '"' <<");\n";
    out_file <<
        "  hist2->SetYTitle(" << '"' << "d(exposure)/d(cos(theta))" << '"' <<");\n";
    
    
    out_file <<
        "  fp = fopen(" << '"' << "./data2.dat" << '"' << "," << '"' << "r" << '"' << ");\n"
        
        "  while ( fscanf(fp," << '"' << "%f" << '"' << ",&ptmp) != EOF ){\n"
        "    p[i++]=ptmp;\n"
        "    if (i==2){\n"
        "      i=0; \n"
        "      iline++;\n"
        "      hist2->Fill(p[0],p[1]); \n"
        "    }\n"
        "  }\n"
        
        
        "  hist2.Draw("
        //<< '"' << "COLZ" << '"' 
        << ");\n"
        
        "}\n";
    
    out_file.close();
    
    system("root -l graph2.cxx");

}
//------------------------------------------------------------
void FluxTestAlg::displayExposure(){
    //make the file
    std::ofstream out_file("data.dat", std::ios::ate);
    
    int i,j;
    for(i=0 ; i<360 ; i++){
        //std::strstream out;
        for(j=0 ; j<180 ; j++){
            //out << m_exposedArea[i][j] << " ";
            out_file << i << "  " << j-90 << "  " << m_exposedArea[i][j] << std::endl;
        }
        //out << std::endl;
        //std::cout << out.str();
    }
    out_file.close();
    //then use rootDisplay to display the stuff in the file
    rootDisplay();
}
//------------------------------------------------------------------    
void FluxTestAlg::rootDisplay(){
    
    std::ofstream out_file("graph.cxx", std::ios::app);
    
    out_file.clear();
    
    out_file << 
        "{\n"
        
        "  FILE *fp;\n"
        "  float ptmp,p[20];\n"
        "  int i, iline=0;\n";
    
    
    out_file <<
        "  TH2D *hist1 = new TH2D(" << '"' << "hist1" << '"' << "," << '"' << "Total Exposure" << '"' << ",360,0.,360.,180,-90.,90.);\n";
    
    out_file <<
        "  hist1->SetXTitle(" << '"' << "l" << '"' <<");\n";
    out_file <<
        "  hist1->SetYTitle(" << '"' << "b" << '"' <<");\n";
    
    
    out_file <<
        "  fp = fopen(" << '"' << "./data.dat" << '"' << "," << '"' << "r" << '"' << ");\n"
        
        "  while ( fscanf(fp," << '"' << "%f" << '"' << ",&ptmp) != EOF ){\n"
        "    p[i++]=ptmp;\n"
        "    if (i==3){\n"
        "      i=0; \n"
        "      iline++;\n"
        "      hist1->Fill(p[0],p[1],p[2]); \n"
        "    }\n"
        "  }\n"
        
        
        "  hist1.Draw("
        << '"' << "COLZ" << '"' 
        << ");\n"
        
        "}\n";
    
    out_file.close();
    
    system("root -l graph.cxx");
    
}




/*
*  This routine will convert longitude and latitude values ("l" and "b")
*  into equivalent cartesian coordinates "x" and "y".  The return values
*  will be in internal Aitoff units where "*x" is in the range [-2, 2]
*  and "*y" is in the range [-1, 1].
*/
std::pair<double,double> FluxTestAlg::hammerAitoff(double l,double b){
    
    double x, y;
    float lover2, den;
    float radl, radb;
    
    /*
    *  Force "l" to be in the range (-180 <= l <= +180) and
    *  "b" to be in the range (-90 <= b <= +90).
    */
    
    while (l < -180.0) l += 360;
    while (l >  180.0) l -= 360;
    while (b <  -90.0) b += 180;
    while (b >   90.0) b -= 180;
    
    /*  Convert l and b to radians. */
    
    double RadPerDeg = M_PI/180.;
    radl = l * RadPerDeg;
    radb = b * RadPerDeg;
    
    lover2 = radl / 2.0;
    den = sqrt(1.0 + (cos(radb) * cos(lover2)));
    x = 2.0 * cos(radb) * sin(lover2) / den;
    y = sin(radb) / den;
    
    /* "x" is now in the range [-2, 2] and "y" in the range [-1, 1]. */
    
    x +=2;
    y +=1;
    x*=90;
    y*=90;
    return std::make_pair<double,double>(x,y);
}
//--------------------------------------------------------------------------------------------
double FluxTestAlg::zenithExposure(double zenithAngle){
    double characteristicAngle=25;
    if(zenithAngle <=60)return pow(2.71828,-zenithAngle/characteristicAngle);
    return 0;
}
double FluxTestAlg::energyExposure(double energy){
    return 1;
}
double FluxTestAlg::azimuthExposure(double azimuth){
    //note:  azimuth here is in the range (-180,180), where 0 represents the x-axis of the spacecraft.
    //if(azimuth<=0){
    return 1.;
    //}else{return 0;}
}

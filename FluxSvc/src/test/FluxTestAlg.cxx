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
    double m_glastExposureAngle;
    int m_exposureMode;
    double m_exposedArea[360][180];
    double m_currentTime;
    double m_passedTime;  //time passed during this event
    std::vector<exposureSet> findExposed(double l,double b, double deltat);
    void addToTotalExposure(std::vector<exposureSet>);
    void displayExposure();
    void rootDisplay();

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
    declareProperty("exposureMode", m_exposureMode=1.);
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
    //fsvc->pass(1000000.);
    
    
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
    MsgStream log(msgSvc(), name());
    
    std::vector<exposureSet> returned;
    double angularRadius = m_glastExposureAngle/2.;
    
    if(m_exposureMode == 0){
        exposureSet point;
        point.x = l; //yes, this is doing an implicit cast.
        point.y = b;
        point.amount = deltat;
        returned.push_back(point);
    }
    
    else if(m_exposureMode == 1){
        for(int i= l-angularRadius ; i<=l+angularRadius ; i++){
            for(int j= b-angularRadius ; j<=b+angularRadius ; j++){
                
                if((pow(l-i,2)+pow(b-j,2) <= pow(angularRadius,2))){
                    //set up the point, and stick it into the vector
                    exposureSet point;// = new exposureSet;
//                    float correctedl = fmod(i, 360);   // Fold into the range [0, 360)
//                    float correctedb = fmod(j, 180);   // Fold into the range [0, 360)
//                    if(correctedl < 0)correctedl+=360;
//                    if(correctedb < 0)correctedb+=180;
                      float correctedl = i;   // Fold into the range [0, 360)
                      float correctedb = j;   // Fold into the range [0, 360)
                   
                      if(correctedl < 0)correctedl+=360;
                      if(correctedb < /*-90*/0)correctedb+=180;
                      if(correctedl > 360)correctedl-=360;
                      if(correctedb > /*90*/180)correctedb-=180;


                    point.x = correctedl; //yes, this is doing an implicit cast.
                    point.y = correctedb;
                    point.amount = deltat;
                    returned.push_back(point);
                }
            }
        }
    }
    else if(m_exposureMode == 2){
      for(int i= l-angularRadius ; i<=l+angularRadius ; i++){
            for(int j= b-angularRadius ; j<=b+angularRadius ; j++){
                
                if((pow(l-i,2)+pow(b-j,2) <= pow(angularRadius,2))){
                    //set up the point, and stick it into the vector
                    exposureSet point;
                    float correctedl = fmod(i, 360);   // Fold into the range [0, 360)
                    float correctedb = fmod(j, 180);   // Fold into the range [0, 360)
                    if(correctedl < 0)correctedl+=360;
                    if(correctedb < 0)correctedb+=180;

                    point.x = correctedl; //yes, this is doing an implicit cast.
                    point.y = correctedb;

                    //here, the amount depends on the time, and a gaussian in azimuth
                    double gaussian = pow(2.71828, -(pow(l-i,2)+pow(b-j,2))/(2*pow(angularRadius,2)));
                    point.amount = deltat*gaussian;
                    returned.push_back(point);
                }
            }
        }
    }
    else if(m_exposureMode == 3){
        for(int i= 0 ; i<360 ; i++){
            for(int j= 0 ; j<180 ; j++){
                double degToRad = M_PI/180;
                
                if(1.5*acos(sin((j-90.)*degToRad)*sin((b-90.)*degToRad)+cos((b-90.)*degToRad)*cos((j-90.)*degToRad)*cos((l-i)*degToRad))/degToRad  <= pow(angularRadius,1)){
                    //set up the point, and stick it into the vector
                    exposureSet point;// = new exposureSet;
                      float correctedl = i;   // Fold into the range [0, 360)
                      float correctedb = j;   // Fold into the range [0, 360)
                   
                      if(correctedl < 0)correctedl+=360;
                      if(correctedb < /*-90*/0)correctedb+=180;
                      if(correctedl > 360)correctedl-=360;
                      if(correctedb > /*90*/180)correctedb-=180;


                    point.x = correctedl; //yes, this is doing an implicit cast.
                    point.y = correctedb;
                    point.amount = deltat;
                    returned.push_back(point);
                }
            }
        }
    }
        else if(m_exposureMode == 4){
        for(int i= l-angularRadius ; i<=l+angularRadius ; i++){
            for(int j= b-angularRadius ; j<=b+angularRadius ; j++){
                
                if((pow(l-i,2)+pow(b-j,2) <= pow(angularRadius,2))){
                    //set up the point, and stick it into the vector
                    exposureSet point;// = new exposureSet;
//                    float correctedl = fmod(i, 360);   // Fold into the range [0, 360)
//                    float correctedb = fmod(j, 180);   // Fold into the range [0, 360)
//                    if(correctedl < 0)correctedl+=360;
//                    if(correctedb < 0)correctedb+=180;
                      float correctedl = i;   // Fold into the range [0, 360)
                      float correctedb = j;   // Fold into the range [0, 360)
                   
                      if(correctedl < 0)correctedl+=360;
                      if(correctedb < /*-90*/0)correctedb+=180;
                      if(correctedl > 360)correctedl-=360;
                      if(correctedb > /*90*/180)correctedb-=180;

                      std::pair<double,double> abc = hammerAitoff(correctedl,correctedb);
                    point.x = abc.first; //yes, this is doing an implicit cast.
                    point.y = abc.second;
                    point.amount = deltat;
                    returned.push_back(point);
                }
            }
        }
    }

    else{
        log << MSG::ERROR << "Invalid Exposure Mode " << endreq;
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
            if(x >=0 && y>=0 && x<360 && y<180) m_exposedArea[x][y] += amount;
            
        }
    }else{
        
        std::cout << "error in addToTotalExposure - null vector input" <<std::endl;
    }
}

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

void FluxTestAlg::rootDisplay(){
    
    std::ofstream out_file("graph.cxx", std::ios::app);
    
    out_file.clear();
    
    out_file << 
        "{\n"
        
        "  FILE *fp;\n"
        "  float ptmp,p[20];\n"
        "  int i, iline=0;\n";
    out_file << 
        "  ntuple = new TNtuple(" << '"' << "ntuple" << '"' << "," << '"' << "NTUPLE" << '"' << "," << '"' <<"x:y:z" << '"' << ");\n";
    
    
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
        "      ntuple->Fill(p[0],p[1],p[2]); \n"
        "      hist1->Fill(p[0],p[1],p[2]); \n"
        "    }\n"
        "  }\n"
        
        
        //"  ntuple.Draw(" << '"' << "z:y:x" << '"' 
        //<< "," << '"' << '"'
        //<< "," << '"' << "HIST" << '"' 
        //<< ");\n"
        
        "  hist1.Draw("
        << '"' << "COLZ" << '"' 
        << ");\n"
        
        "}\n";
    
    out_file.close();
    
    system("root -l graph.cxx");
    
}




/////////T$EST FILE////////////////////////////
/*
 *  This routine will convert longitude and latitude values ("l" and "b")
 *  into equivalent cartesian coordinates "x" and "y".  The return values
 *  will be in internal Aitoff units where "*x" is in the range [-2, 2]
 *  and "*y" is in the range [-1, 1].
 */
//#ifdef PROTOTYPE
//static void aitoffConvert(float l, float b, float *x, float *y)
//#else
std::pair<double,double> FluxTestAlg::hammerAitoff(double l,double b){

//#endif /* PROTOTYPE */
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

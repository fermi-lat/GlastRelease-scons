/** @file FesExposure.cxx
    @brief declare and implement the Algorithm FesExposure

    $Header$

*/
// Include files

// Gaudi system includes
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"

#include "Event/MonteCarlo/Exposure.h"

#include "astro/SkyDir.h"
#include "astro/EarthCoordinate.h"
#include "astro/SolarSystem.h"

//flux
#include "FluxSvc/IFluxSvc.h"
#include "flux/IFlux.h"
#include "astro/GPS.h"


#include <cassert>
#include <vector>
#include <fstream>
#include <iomanip>

// Record Size LDF Header (4 32bit) +
//             5 Attitude 5*(13 32bit) + 
//             1 Position (8 32 bit) = 77 32bit words 
#define attRecordSize 77
// Time from Jan 1, 2001 to launch  ~7 years * 30 million seconds
//THB time is already in MET #define launchOffset  210000000
/** 
* \class FesExposure
*
* \brief This is an Algorithm designed to get information about the GLAST locatation and orientaion
* and save it in a special tuple
*
* \author Brian Winer, mods by T. Burnett
* 
*/
    static inline unsigned int writeAtt (unsigned int* data, FILE* fp);
    static inline void swap (unsigned int *wrds, int nwrds);
    static inline Hep3Vector fromRaDec(double RA, double Dec);

class FesExposure : public Algorithm {
public:
    FesExposure(const std::string& name, ISvcLocator* pSvcLocator);

    //stuff that an Algorithm needs.
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    StatusCode addAttitudeHeader();
    StatusCode calcAttitudeContribution(CLHEP::Hep3Vector xaxis, CLHEP::Hep3Vector zaxis);
    StatusCode calcPositionContribution(double x, double y, double z);
    StatusCode convertTime(double time);

private: 

    int m_tickCount;
    IFluxSvc*   m_fluxSvc;

    unsigned int* m_att_data_start;
    unsigned int* m_att_data;
    int           m_attCont;
    FILE*         m_outfile;
    bool          m_startCont;
    unsigned int  m_secTime;
    unsigned int  m_microTime;
    double        m_lastx,m_lasty,m_lastz;
    double        m_lastAttTime, m_lastPosTime;
    Hep3Vector    m_lastXaxis,m_lastYaxis,m_lastZaxis;
    std::string   m_FileName;
    StringProperty m_timerName;
    
    typedef union IntDble_t {
       unsigned int dui[2];
       double dbl;
    } IntDble;

    typedef union IntFlt_t {
       unsigned int intVal;
       float fltVal;
    } IntFlt;
    
    
};
//------------------------------------------------------------------------

static const AlgFactory<FesExposure>  Factory;
const IAlgFactory& FesExposureFactory = Factory;

//------------------------------------------------------------------------
//! ctor
FesExposure::FesExposure(const std::string& name, ISvcLocator* pSvcLocator)
: Algorithm(name, pSvcLocator)
, m_tickCount(0)
, m_attCont(0)
, m_outfile(0)
, m_startCont(true)
, m_secTime(0)
, m_microTime(0)
{
  declareProperty("FileName"  ,m_FileName="Attitude");
  declareProperty("TimerName" ,m_timerName="timer_5Hz");

}

//------------------------------------------------------------------------
//! set parameters and attach to various perhaps useful services.
StatusCode FesExposure::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    m_att_data_start = 0;

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

   /* Allocate a buffer big enough to hold a maximally sized event */
   unsigned char* ptr = (unsigned char *)malloc (attRecordSize*4);
   if (ptr == 0) { return StatusCode::FAILURE;  }
   m_att_data_start = (unsigned int *)ptr;
   m_att_data = m_att_data_start;
    

// Open up the output file
   char file[120];
   sprintf(file,"%s.ldf",m_FileName.c_str());

   m_outfile = fopen (file, "wb");
   if (!m_outfile)
   {
     /* Error in opening the output file.. */
     free (ptr);
     log << MSG::ERROR << "Unable to open file to write Attitude Information" << endreq;
     return StatusCode::FAILURE;
   }    

    m_lastx = m_lasty = m_lastz = 0.;
    m_lastAttTime = m_lastPosTime = 0.;
    m_lastXaxis = Hep3Vector(0.,0.,0.);
    m_lastYaxis = Hep3Vector(0.,0.,0.);
    m_lastZaxis = Hep3Vector(0.,0.,0.);
   
    return sc;
}

static void  swap (unsigned int *wrds, int nwrds)
/*
   DESCRIPTION
   -----------
   Performs a byte-swap on the specified number of 32-bit words. This
   is an El-Cheapo implementation, it should not be used for heavy
   duty byte-swapping.

   PARAMETERS
   ----------
         wrds: The 32-bit words to byte-swap

        nwrds: The number of 32-bit words to byte-swap. 

   RETURNS
   -------
   Nothing
*/
{
   while (--nwrds >= 0)
   {
     unsigned int tmp;
     
     tmp = *wrds;
     tmp = ((tmp & 0x000000ff) << 24) |
       ((tmp & 0x0000ff00) <<  8) |
       ((tmp & 0x00ff0000) >>  8) |
       ((tmp & 0xff000000) >> 24); 
     *wrds++ = tmp;
   }
   
   return;
}

static unsigned int writeAtt (unsigned int* data, FILE* fp)
{
   /*
    |  !!! KLUDGE !!!
    |  --------------
    |  Need to know whether the machine executing this code is a big
    |  or little endian machine. I've checked around and found no one
    |  who can tell me of a compiler defined symbol containing this
    |  tidbit of information. Currently, I've only every run GLEAM on
    |  Intel processors, so I've hardwired this to be little-endian.
   */
    swap (data, attRecordSize);
    if(fp) return fwrite (data, sizeof (*data), attRecordSize, fp);
    return 0;
}

//------------------------------------------------------------------------
//! process an event
StatusCode FesExposure::execute()
{
    using astro::GPS;
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    IFlux* flux=m_fluxSvc->currentFlux();
        
    if(flux->name() != m_timerName.value()) return sc;


    //printf("found time tick\n");
    // The GPS singleton has current time and orientation
    GPS* gps = m_fluxSvc->GPSinstance();
    double time = gps->time();
    CLHEP::Hep3Vector location = 1.e3* gps->position(); // special, needs its own time
    
    // cartesian location of the LAT (in m)
    double posX = location.x();
    double posY = location.y();
    double posZ = location.z(); 

/*  printf("FES::Attitude Information:  
             Time:     %f
             Latitude: %f
             Longitude:%f
             Altitude: %f
             PosX:     %f
             PosY:     %f
             PosZ:     %f
             RightAsX: %f
             DelX:     %f
             RightAsZ: %f
             DelZ:     %f \n",time,latitude,longitude,altitude,posX,posY,posZ,RightAsX,DeclX,RightAsZ,DeclZ);      
*/

// Convert Time to seconds and microseconds since Jan 1, 2001
    convertTime(time);

// If this is a new contribution, add the header
   if(m_startCont) {
      m_att_data = m_att_data_start;
      addAttitudeHeader();
//      printf("Adding the Header\n");
      m_startCont = false;
   }

// Calculate the attitude Contribution
   calcAttitudeContribution(gps->xAxisDir()(), gps->zAxisDir()());
//   printf("Adding the Att. Contribution # %i at time %f\n",m_attCont,time);
   m_attCont++;

// If it is the end of 1 second, calculate the position information
   if(m_attCont > 4) { 
     calcPositionContribution(posX,posY,posZ);
//     printf("Adding the Position Cont at time %f\n",time);
   
// flush this contributions to the file
     writeAtt(m_att_data_start,m_outfile);
//     printf("Dumping the Att file\n");
     
// Next time through we need to start a new contribution     
     m_startCont = true;
     m_attCont = 0;
   }

    m_tickCount++;
    return sc;
}

StatusCode FesExposure::addAttitudeHeader() {

   StatusCode sc = StatusCode::SUCCESS;
   
   int LATdatagramID = 0xF0002;
   int LD_Version    = 0x800; //MSB must be set
   *m_att_data++ = (LD_Version)<<20 | (LATdatagramID & 0xfffff);    
   *m_att_data++ = (attRecordSize & 0xffffff);
       
   int LATcontribID = 0xFFFFF;
   int LC_Version   = 0x000;
   *m_att_data++ = (LC_Version)<<20 | (LATcontribID & 0xfffff);
   *m_att_data++ = ((attRecordSize-2) & 0xffffff); 

   
   return sc;
}

StatusCode FesExposure::calcAttitudeContribution(Hep3Vector xAxis, Hep3Vector zAxis){
   StatusCode sc = StatusCode::SUCCESS;

// Now get the directions starting from the Z axis. Why?
   Hep3Vector yAxis = zAxis.cross(xAxis).unit();
   xAxis = yAxis.cross(zAxis);

// Now Convert these directional vectors to a Quaternion
   double s = 1.0 + xAxis.x() + yAxis.y() + zAxis.z();
   IntDble q1,q2,q3,q4;
   if(s > 0.01) {
      q1.dbl = yAxis.z() - zAxis.y();
      q2.dbl = zAxis.x() - xAxis.z();
      q3.dbl = xAxis.y() - yAxis.x();
      q4.dbl = s;   
   } else {
      if(xAxis.x() >= yAxis.y() && xAxis.x() >= zAxis.z()) {
        q1.dbl = 1.0 + xAxis.x() - yAxis.y() - zAxis.z();
        q2.dbl = yAxis.x() + xAxis.y();
        q3.dbl = zAxis.x() + xAxis.z();
        q4.dbl = yAxis.z() - zAxis.y();  
      } else if(yAxis.y() >= xAxis.x() && yAxis.y() >= zAxis.z()) {
        q1.dbl = yAxis.x() + xAxis.y();
        q2.dbl = 1.0 - xAxis.x() + yAxis.y() - zAxis.z();
        q3.dbl = zAxis.y() + yAxis.z();
        q4.dbl = zAxis.x() - xAxis.z();  
      } else {
        q1.dbl = zAxis.x() + xAxis.z();
        q2.dbl = zAxis.y() + yAxis.z();
        q3.dbl = 1.0 - xAxis.x() - yAxis.y() + zAxis.z();
        q4.dbl = xAxis.y() - yAxis.x();  
      }
   
   }    

// Determine the normalization of the components and then scale so
// the quaternion has unit length.
   double norm = q1.dbl*q1.dbl+q2.dbl*q2.dbl+q3.dbl*q3.dbl+q4.dbl*q4.dbl;
   if(norm >0.) {
      norm = sqrt(norm);
      q1.dbl = q1.dbl/norm;
      q2.dbl = q2.dbl/norm;
      q3.dbl = q3.dbl/norm;
      q4.dbl = q4.dbl/norm;
   }
/*
   printf("RaX %f  DecX %f  RaX %f  DecX %f\n",RAx,DecX,RAz,DecZ);
   printf("X-Axis Vector i,j,k,mag: %f %f %f %f Angle y %f Angle z %f\n",xAxis.x(),xAxis.y(),xAxis.z(),xAxis.mag(),xAxis.angle(yAxis),xAxis.angle(zAxis));
   printf("Y-Axis Vector i,j,k,mag: %f %f %f %f\n",yAxis.x(),yAxis.y(),yAxis.z(),yAxis.mag());
   printf("Z-Axis Vector i,j,k,mag: %f %f %f %f\n",zAxis.x(),zAxis.y(),zAxis.z(),zAxis.mag());
   printf("Q1 %f Q2 %f Q3 %f Q4 %f\n",q1.dbl,q2.dbl,q3.dbl,q4.dbl);
*/   
   // First put in the time stamp
   *m_att_data++ = m_secTime;
   *m_att_data++ = m_microTime;
   
   // Now Quaternion q1 (64 bits)
   *m_att_data++ = q1.dui[1];
   *m_att_data++ = q1.dui[0];

   // Now Quaternion q2 (64 bits)
   *m_att_data++ = q2.dui[1];
   *m_att_data++ = q2.dui[0];

   // Now Quaternion q3 (64 bits)
   *m_att_data++ = q3.dui[1];
   *m_att_data++ = q3.dui[0];

   // Now Quaternion q4 (64 bits)
   *m_att_data++ = q4.dui[1];
   *m_att_data++ = q4.dui[0];

   
   // Now Angular Velocity x,y,z
   // This is not part of Gleam/GPS(?) so extrapolate from
   // previous point.  This info is not used in FSW
   // First entry fill with zeros
   // Time since last Attitude Info 
     double newTime = (double)m_secTime + ((double)m_microTime)/1000000.;
     double deltaTime = newTime - m_lastAttTime;

   if(m_lastXaxis.mag() == 0) { 

       IntFlt dummy;
       dummy.fltVal = 0.;
       *m_att_data++ = dummy.intVal; 
       *m_att_data++ = dummy.intVal; 
       *m_att_data++ = dummy.intVal;

   } else {
       
// get change in direction of axes
     IntFlt angVelX,angVelY,angVelZ;
     if( deltaTime>0) {
         angVelX.fltVal = xAxis.angle(m_lastXaxis)/deltaTime;
         angVelY.fltVal = yAxis.angle(m_lastYaxis)/deltaTime;
         angVelZ.fltVal = zAxis.angle(m_lastZaxis)/deltaTime;   
     }
     else{
        angVelX.intVal=angVelY.intVal=angVelZ.intVal = 0;
     }
//     printf("Angular Velocity (rad/sec) x,y,z,dt: %f %f %f %f\n",angVelX.fltVal,angVelY.fltVal,angVelZ.fltVal,deltaTime);
     
// Stuff the Angular Velocity into the Record   
     *m_att_data++ = angVelX.intVal; 
     *m_att_data++ = angVelY.intVal; 
     *m_att_data++ = angVelZ.intVal;
     
   }     
   
     m_lastXaxis = xAxis;
     m_lastYaxis = yAxis;
     m_lastZaxis = zAxis;
     m_lastAttTime = newTime;

   
   return sc;
}

static inline Hep3Vector fromRaDec(double RA, double Dec){

    double theta    = (90.0 - Dec)*M_PI/180.0;
    double phi      = RA*M_PI/180.0;
    double sinTheta = sin(theta);
       
    Hep3Vector vec = Hep3Vector(sinTheta*cos(phi), sinTheta*sin(phi), cos(theta));
 
    return vec;
}

StatusCode FesExposure::calcPositionContribution(double x, double y, double z) {

   StatusCode sc = StatusCode::SUCCESS;

   // First put in the time stamp
   *m_att_data++ = m_secTime;
   *m_att_data++ = m_microTime;

   // Setup variables
   IntFlt xpos,ypos,zpos;
   xpos.fltVal = (float)x;
   ypos.fltVal = (float)y;
   zpos.fltVal = (float)z;

   // Now Position x,y,z
   *m_att_data++ = xpos.intVal; 
   *m_att_data++ = ypos.intVal; 
   *m_att_data++ = zpos.intVal;      

   // Now Velocity  x,y,z
   // this info is not really in GLEAM/GPS(?) For now make an
   // estimate based on the last postion.  The Flight software 
   // does really make use of this information
   double newTime = (double)m_secTime + ((double)m_microTime)/1000000.;
   double deltaTime = newTime - m_lastPosTime;
   
   IntFlt vx,vy,vz;
   if(m_lastx == 0. && m_lasty == 0. && m_lastz == 0.) {
      vx.fltVal = vy.fltVal = vz.fltVal = 0.;   
   } else {   
      vx.fltVal = (float)(x-m_lastx)/deltaTime;
      vy.fltVal = (float)(y-m_lasty)/deltaTime;
      vz.fltVal = (float)(z-m_lastz)/deltaTime;
   }
   
   // Putting in the estimate of the velocity
   *m_att_data++ = vx.intVal; 
   *m_att_data++ = vy.intVal; 
   *m_att_data++ = vz.intVal;      
   
   // Store for next time through
   m_lastPosTime = newTime;
   m_lastx = x;
   m_lasty = y;
   m_lastz = z;
   
//   printf("Position x,y,z  %f, %f, %f\n",x,y,z);
//   printf("Velocity x,y,z  %f, %f, %f\n",vx.fltVal,vy.fltVal,vz.fltVal);
   
   return sc;
}

//------------------------------------------------------------------------
//! clean up, summarize
StatusCode FesExposure::convertTime(double time){
    // finish up

    StatusCode  sc = StatusCode::SUCCESS;

    static int launchOffset(0);  // THB since time is already in MET    
    m_secTime = launchOffset + (int)time;
    m_microTime = (int)(1000000.*(time - (int)time));
//    printf("Time %f seconds %i microseconds %i\n",time,m_secTime,m_microTime);

    return sc;
}


//------------------------------------------------------------------------
//! clean up, summarize
StatusCode FesExposure::finalize(){
    // finish up

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "FesExposure Processed " << m_tickCount << " ticks" << endreq;

    if(m_outfile) fclose(m_outfile);

    return sc;
}


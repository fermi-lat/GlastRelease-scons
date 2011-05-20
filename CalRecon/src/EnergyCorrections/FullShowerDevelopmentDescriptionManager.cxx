/**   
*
* THIS CODE IS STILL NOT FULLY DOCUMENTED !!!!!!!!!!!!!!!!!!!!
* DESCRIPTION COMMENTS WILL BE ADDED SOON.
* Philippe Bruel
*
*/

#include "FullShowerDevelopmentDescriptionManager.h"

#include "idents/TowerId.h"
#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Transform3D.h"


#include "TMath.h"

FullShowerGeometryManager::FullShowerGeometryManager(IGlastDetSvc *m_detSvc_input)
  :m_detSvc(m_detSvc_input)
{
  Initialize();
}

FullShowerGeometryManager::~FullShowerGeometryManager()
{
  
}

void FullShowerGeometryManager::Initialize()
{
  // Retrieve geometry from Detector Service

  FSGM_CsIHeight = 19.9;
  m_detSvc->getNumericConstByName(std::string("CsIHeight"),&FSGM_CsIHeight);
  
  FSGM_CsIWidth = 26.7;
  m_detSvc->getNumericConstByName(std::string("CsIWidth"),&FSGM_CsIWidth);
  
  FSGM_CsILength = 326.0;
  m_detSvc->getNumericConstByName(std::string("CsILength"),&FSGM_CsILength);

  FSGM_cellVertPitch = 21.35;
  m_detSvc->getNumericConstByName(std::string("cellVertPitch"),&FSGM_cellVertPitch);

  FSGM_cellHorPitch = 27.84;
  m_detSvc->getNumericConstByName(std::string("cellHorPitch"),&FSGM_cellHorPitch);

  FSGM_towerPitch = 374.5;
  m_detSvc->getNumericConstByName(std::string("towerPitch"),&FSGM_towerPitch);

  // FSGM_crackHalfWidth is the average of X and Y layers crackHalfWidths
  FSGM_crackHalfWidth = ((FSGM_towerPitch-FSGM_CsILength)/2 + (FSGM_towerPitch-12*FSGM_cellHorPitch)/2)/2;

  // Retrieving the top position of the calorimeter as in TkrUtil/*/src/TkrGeometrySvc.cxx
  FSGM_calZTop = -47.395;

  int m_Tower = 0;
  idents::VolumeIdentifier topLayerId;
  topLayerId.init(0,0);
  topLayerId.append(0);               // in Tower
  idents::TowerId t(m_Tower);  
  topLayerId.append(t.iy());          // yTower
  topLayerId.append(t.ix());          // xTower
  topLayerId.append(0);  // CAL
  topLayerId.append(0);  // layer
  topLayerId.append(0);  // x view
  StatusCode sc;
  HepGeom::Transform3D transfTop;
  int count;
  for (count=0;count<3;++count) {
    topLayerId.append(0);
    if((sc = m_detSvc->getTransform3DByID(topLayerId,&transfTop)).isSuccess()) break;
  }
  if(sc.isSuccess()) 
    {
      FSGM_calZTop = (transfTop.getTranslation()).z()+ FSGM_cellVertPitch/2;
    }

  FSGM_calZBot = FSGM_calZTop-8*FSGM_cellVertPitch;

  // prepare safe volume in order to define the segment of the trajectory used by the propagator
  double safedistance = 50;
  FSGM_CalSafeBoundaries[0][0] = -2*FSGM_towerPitch-safedistance;
  FSGM_CalSafeBoundaries[0][1] = +2*FSGM_towerPitch+safedistance;
  FSGM_CalSafeBoundaries[1][0] = -2*FSGM_towerPitch-safedistance;
  FSGM_CalSafeBoundaries[1][1] = +2*FSGM_towerPitch+safedistance;
  FSGM_CalSafeBoundaries[2][0] = FSGM_calZBot-safedistance;
  FSGM_CalSafeBoundaries[2][1] = FSGM_calZTop+safedistance;

  // get the number of towers
  m_detSvc->getNumericConstByName("xNum", &FSGM_numX);
  m_detSvc->getNumericConstByName("yNum", &FSGM_numY);
  FSGM_flight_geom = true;
  if(FSGM_numX==4 && FSGM_numY==1) FSGM_flight_geom = false;

  // Fill the look-up table for even layers
  geom_xy_step = (4*FSGM_towerPitch)/(double)(FSGM_XY_MAX);
  
  double xyz[3];
  int whereincal[4];

  int i,j;
  // reference = middle of first layer
  xyz[2] = FSGM_calZTop-FSGM_cellVertPitch/2;
  
  // reference = middle of a log
  xyz[1] = FSGM_towerPitch/2+FSGM_cellHorPitch/2;
  for(i=0;i<FSGM_XY_MAX;++i)
    {
      xyz[0] = -2*FSGM_towerPitch + geom_xy_step*(0.5+(double)i);
      WhereInCalForGeom(xyz,whereincal);
      switch(whereincal[0])
        {
        case 0:
          geom_x_mat[i][0] = false;
          geom_x_mat[i][1] = false;
          break;
        case 1:
          geom_x_mat[i][0] = true;
          geom_x_mat[i][1] = false;
          break;
        case 2:
          geom_x_mat[i][0] = false;
          geom_x_mat[i][1] = true;
          break;
        default:
          geom_x_mat[i][0] = false;
          geom_x_mat[i][1] = false;
          break;
        }
    }
  
  // reference = middle of a log
  xyz[0] = FSGM_towerPitch/2+ FSGM_cellHorPitch/2;
  for(j=0;j<FSGM_XY_MAX;++j)
    {
      xyz[1] = -2*FSGM_towerPitch+geom_xy_step*(0.5+(double)j);
      WhereInCalForGeom(xyz,whereincal);
      switch(whereincal[0])
        {
        case 0:
          geom_y_mat[j][0] = false;
          geom_y_mat[j][1] = false;
          break;
        case 1:
          geom_y_mat[j][0] = true;
          geom_y_mat[j][1] = false;
          break;
        case 2:
          geom_y_mat[j][0] = false;
          geom_y_mat[j][1] = true;
          break;
        default:
          geom_y_mat[j][0] = false;
          geom_y_mat[j][1] = false;
          break;
        }
    }

  // Fill the test points used by the propagator
  FillXYPoints();

  // Parameters used for the radial shower development description
  FSGM_CORE0 = 0.0211192;
  FSGM_CORE1 = 0.114962;
  FSGM_TAIL0 = 0.450839;
  FSGM_TAIL1 = 0.840738;
  FSGM_TAIL2 = -1.87526;
  FSGM_TAIL3 = 0.643606;
  FSGM_PCT0 = 2.51719;
  FSGM_PCT1 = 0.554606;
  FSGM_PCT2 = 0.851856;

  // Limits for the relative position in the shower to shower maximum
  FSGM_tmin = 0.25;
  FSGM_tmax = 2.5;

  // Fill radial profile look-up table
  FillRadialProfile();

}

void FullShowerGeometryManager::FillXYPoints()
{
  int n = (int)sqrt(1/3.2*(double)FSGM_NPOINTS_MAX);
  int i,j;
  double step = 1/(double)(n);
  double x,y;
  FSGM_NCircle = 0;
  for(i=0;i<=2*n;++i)
    {
      x = -1+step*(double)i;
      for(j=0;j<=2*n;++j)
        {
          y = -1+step*(double)j;
          if(sqrt(x*x+y*y)>1) continue;
          FSGM_XCircle[FSGM_NCircle] = x;
          FSGM_YCircle[FSGM_NCircle] = y;
          FSGM_RCircle[FSGM_NCircle] = sqrt(x*x+y*y);
          ++FSGM_NCircle;
        }
    }
}

void FullShowerGeometryManager::WhereInCalForGeom(double *xyz, int *whereincal)
{
  whereincal[0] = 0; // material
  whereincal[1] = 0; // layer
  whereincal[2] = 0; // tower
  whereincal[3] = 0; // column
  if(xyz[2]>=FSGM_calZTop) return;
  if(xyz[2]<FSGM_calZBot) return;
  int ilayer = (int)floor( (FSGM_calZTop-xyz[2])/FSGM_cellVertPitch );
  //  if(fabs(xyz[2]- (FSGM_calZTop-FSGM_cellVertPitch*(0.5+(double)ilayer)))>FSGM_CsIHeight*0.5) return;
  //
  if(fabs(xyz[0])>2*FSGM_towerPitch) return;
  if(fabs(xyz[1])>2*FSGM_towerPitch) return;
  whereincal[2] = (int)floor((xyz[0]+2*FSGM_towerPitch)/FSGM_towerPitch);
  whereincal[3] = (int)floor((xyz[1]+2*FSGM_towerPitch)/FSGM_towerPitch);
  whereincal[2] = 4*whereincal[3]+whereincal[2];
  //
  double xmod = xyz[0]+2*FSGM_towerPitch-FSGM_towerPitch*floor((xyz[0]+2*FSGM_towerPitch)/FSGM_towerPitch);
  double ymod = xyz[1]+2*FSGM_towerPitch-FSGM_towerPitch*floor((xyz[1]+2*FSGM_towerPitch)/FSGM_towerPitch);
  if(ilayer%2==0)
    whereincal[3] = (int)floor((ymod-(FSGM_towerPitch/2-6*FSGM_cellHorPitch))/FSGM_cellHorPitch);
  else
    whereincal[3] = (int)floor((xmod-(FSGM_towerPitch/2-6*FSGM_cellHorPitch))/FSGM_cellHorPitch);
  //
  if(FSGM_towerPitch-xmod<xmod) xmod = FSGM_towerPitch-xmod;
  if(FSGM_towerPitch-ymod<ymod) ymod = FSGM_towerPitch-ymod;
  if(xmod<FSGM_crackHalfWidth || ymod<FSGM_crackHalfWidth)
    {
      whereincal[0] = 2;
      return;
    }
  double xy = ymod;
  if(ilayer%2==1) xy = xmod;
  xy -= (FSGM_towerPitch/2-6*FSGM_cellHorPitch);
  double xymod = xy-FSGM_cellHorPitch*floor(xy/FSGM_cellHorPitch);
  if(FSGM_cellHorPitch-xymod<xymod) xymod = FSGM_cellHorPitch-xymod;
  if(xymod<(FSGM_cellHorPitch-FSGM_CsIWidth)/2) return;
/*   double zlayer = (FSGM_calZTop-xyz[2])-FSGM_cellVertPitch*(double)ilayer; */
/*   if(FSGM_cellVertPitch-zlayer>zlayer) zlayer = FSGM_cellVertPitch-zlayer; */
/*   if(zlayer<(FSGM_cellVertPitch-FSGM_CsIHeight)/2) return; */
  whereincal[0] = 1;
  whereincal[1] = ilayer;
}

void FullShowerGeometryManager::WhereInCalForGeomCU(double *xyz, int *whereincal)
{
  whereincal[0] = 0;
  whereincal[1] = 0;
  whereincal[2] = 0;
  whereincal[3] = 0;
  if(xyz[2]>=FSGM_calZTop) return;
  if(xyz[2]<FSGM_calZBot) return;
  int ilayer = (int)floor( (FSGM_calZTop-xyz[2])/FSGM_cellVertPitch );
  if(xyz[0]>2*FSGM_towerPitch || xyz[0]<-FSGM_towerPitch) return;
  if(fabs(xyz[1])>FSGM_towerPitch/2) return;
  whereincal[2] = (int)floor((xyz[0]+2*FSGM_towerPitch)/FSGM_towerPitch);
  double xmod = xyz[0]+2*FSGM_towerPitch-FSGM_towerPitch*floor((xyz[0]+2*FSGM_towerPitch)/FSGM_towerPitch);
  double ymod = xyz[1]+FSGM_towerPitch/2;
  if(ilayer%2==0)
    whereincal[3] = (int)floor((ymod-(FSGM_towerPitch/2-6*FSGM_cellHorPitch))/FSGM_cellHorPitch);
  else
    whereincal[3] = (int)floor((xmod-(FSGM_towerPitch/2-6*FSGM_cellHorPitch))/FSGM_cellHorPitch);
  if(FSGM_towerPitch-xmod<xmod) xmod = FSGM_towerPitch-xmod;
  if(FSGM_towerPitch-ymod<ymod) ymod = FSGM_towerPitch-ymod;
  if(xmod<FSGM_crackHalfWidth || ymod<FSGM_crackHalfWidth)
    {
      whereincal[0] = 2;
      return;
    }
  double xy = ymod;
  if(ilayer%2==1) xy = xmod;
  xy -= (FSGM_towerPitch/2-6*FSGM_cellHorPitch);
  double xymod = xy-FSGM_cellHorPitch*floor(xy/FSGM_cellHorPitch);
  if(FSGM_cellHorPitch-xymod<xymod) xymod = FSGM_cellHorPitch-xymod;
  if(xymod<(FSGM_cellHorPitch-FSGM_CsIWidth)/2) return;
/*   double zlayer = (FSGM_calZTop-xyz[2])-FSGM_cellVertPitch*(double)ilayer; */
/*   if(FSGM_cellVertPitch-zlayer>zlayer) zlayer = FSGM_cellVertPitch-zlayer; */
/*   if(zlayer<(FSGM_cellVertPitch-FSGM_CsIHeight)/2) return; */
  whereincal[0] = 1;
  whereincal[1] = ilayer;
}

void FullShowerGeometryManager::WhereInCalLT(double *xyz, int *whereincal)
{
  whereincal[0] = 0;
  whereincal[1] = 0;
  // tower and column currently not used so not filled for the LT case
  whereincal[2] = 0; // tower
  whereincal[3] = 0; // column

  // test if outside z_dimension
  if(xyz[2]>=FSGM_calZTop) return;
  if(xyz[2]<FSGM_calZBot) return;

  // test if outside xy_dimensions
  if(fabs(xyz[0])>=2*FSGM_towerPitch) return;
  if(fabs(xyz[1])>=2*FSGM_towerPitch) return;

  // get layer number
  int ilayer = (int)floor( (FSGM_calZTop-xyz[2])/FSGM_cellVertPitch );
  whereincal[1] = ilayer;

  // retrieve position in look-up table
  int i = (int)((xyz[0]+2*FSGM_towerPitch)/geom_xy_step);
  int j = (int)((xyz[1]+2*FSGM_towerPitch)/geom_xy_step);
  
  if(ilayer%2==0)  // look-up table set for even layers, invert x and y otherwise
    {
      if(geom_x_mat[i][1] || geom_y_mat[j][1])
        whereincal[0] = 2;
      else
        whereincal[0] = (int)(geom_x_mat[i][0] && geom_y_mat[j][0]);
    }
  else
    {
      if(geom_y_mat[i][1] || geom_x_mat[j][1])
        whereincal[0] = 2;
      else
        whereincal[0] = (int)(geom_y_mat[i][0] && geom_x_mat[j][0]);
    }
}

void FullShowerGeometryManager::WhereInCal(double *xyz, int *whereincal)
{
  if(FSGM_flight_geom)
    {
      //      WhereInCalForGeom(xyz,whereincal); // not needed for saturation handling
      WhereInCalLT(xyz,whereincal);
    }
  else
    WhereInCalForGeomCU(xyz,whereincal);
}

double FullShowerGeometryManager::RCore(double t)
{
  double tt = t;
  if(tt<FSGM_tmin) tt = FSGM_tmin;
  if(tt>FSGM_tmax) tt = FSGM_tmax;
  return (FSGM_CORE0+FSGM_CORE1*tt);
}

double FullShowerGeometryManager::RTail(double t)
{
  double tt = t;
  if(tt<FSGM_tmin) tt = FSGM_tmin;
  if(tt>FSGM_tmax) tt = FSGM_tmax;
  return (FSGM_TAIL0*(exp(FSGM_TAIL2*(tt-FSGM_TAIL1))+exp(FSGM_TAIL3*(tt-FSGM_TAIL1))));
}

double FullShowerGeometryManager::PCore(double t)
{
  double tt = t;
  if(tt<FSGM_tmin) tt = FSGM_tmin;
  if(tt>FSGM_tmax) tt = FSGM_tmax;
  double arg = FSGM_PCT0*exp((FSGM_PCT1-tt)/FSGM_PCT2-exp((FSGM_PCT1-tt)/FSGM_PCT2));
  if(arg<0) arg = 0;
  if(arg>1) arg = 1;
  return arg;
}

double FullShowerGeometryManager::RadialProfile(double t, double r)
{
  double rcore2 = RCore(t);
  rcore2 = rcore2*rcore2;
  double rtail2 = RTail(t);
  rtail2 = rtail2*rtail2;
  double pcore = PCore(t);
  return pcore*(2*rcore2/(r*r+rcore2)/(r*r+rcore2)) + (1-pcore)*(2*rtail2/(r*r+rtail2)/(r*r+rtail2));
}

double FullShowerGeometryManager::GetEffectiveRadius(double t, double radialcontainedfraction)
{
  double rcore2 = RCore(t);
  rcore2 = rcore2*rcore2;
  double rtail2 = RTail(t);
  rtail2 = rtail2*rtail2;
  double pcore = PCore(t);
  
  double a = (1-radialcontainedfraction)/rcore2/rtail2;
  double b = (1-radialcontainedfraction)*(1/rcore2+1/rtail2) - pcore/rtail2 - (1-pcore)/rcore2;
  double c = -radialcontainedfraction;
  double delta = b*b-4*a*c;
  double r1 = (-b+sqrt(delta))/2/a;
  r1 = sqrt(r1);

  if(r1>RadProf_r_max) r1 = RadProf_r_max;

  return r1;
}

void FullShowerGeometryManager::FillRadialProfile()
{
  // Defining limits and steps
  RadProf_r_max = 3.;
  RadProf_t_max = 2.5;  
  RadProf_r_step = RadProf_r_max/(double)FSGM_RPROF_R_MAX;
  RadProf_t_step = RadProf_t_max/(double)FSGM_RPROF_T_MAX;

  int i,j;
  double r,t;
  for(i=0;i<FSGM_RPROF_T_MAX;++i)
    {
      t = RadProf_t_step*(0.5+(double)i);
      for(j=0;j<FSGM_RPROF_R_MAX;++j)
        {
          r = RadProf_r_step*(0.5+(double)j);
          RadProf[i][j] = RadialProfile(t,r);
        }
    }
}

double FullShowerGeometryManager::GetRadialProfile(double t, double r)
{
  int i_t = (int)(t/RadProf_t_step);
  if(i_t>=FSGM_RPROF_T_MAX)
    i_t = FSGM_RPROF_T_MAX-1;

  int i_r = (int)(r/RadProf_r_step);
  if(i_r>=FSGM_RPROF_R_MAX)
    i_r = FSGM_RPROF_R_MAX-1;
  
  return RadProf[i_t][i_r];
}


/**   
* @class FullShowerDevelopmentDescription
* @author Philippe Bruel
*
* Tool that describes the shower developement in the calorimeter given
* the length in X0 seen in the tracker and the position of the shower maximum
*
* $Header$
*/

FullShowerDevelopmentDescription::FullShowerDevelopmentDescription(FullShowerGeometryManager *fsgm_input, int type_input, double zstep_input, double radialcontainedfraction_input)
  :m_fsgm(fsgm_input)
{
  Initialize();
  Reset();
  Type = type_input;
  ZStep = zstep_input;
  RadialContainedFraction = radialcontainedfraction_input;
  wideningfactor = 1.;
}

FullShowerDevelopmentDescription::~FullShowerDevelopmentDescription()
{

}

void FullShowerDevelopmentDescription::Initialize()
{

  // Radiation lengths
  FSDD_MOLRAD = 35.;
  FSDD_gCSI = 4.51;
  FSDD_gCRK = 0.607287449392712508; // 10.*2.7/(2.*22.23) // 10mm
  FSDD_XCSI = 18.5*(m_fsgm->FSGM_cellVertPitch/m_fsgm->FSGM_CsIHeight);
  FSDD_XCRK = 395.364666666666722; // 24.01/FSDD_gCRK*10. // *10 to get it in mm

  Reset();
}

void FullShowerDevelopmentDescription::Reset()
{

  int i,j,k;
  NStep = 0;
  ZStep = 0;
  x0maxshower = 0;
  crackmaxfrac = 0;
  x0crackmaxfrac = 0;
  startx0 = 0;
  lastx0 = 0;
  totx0cal = 0;
  totx0crack = 0;
  for(i=0;i<8;++i) totx0lay[i] = 0;
  for(i=0;i<8;++i) posx0lay[i] = 0;
  for(i=0;i<FSDD_NSTEPS_MAX;++i)
    {
      dX0[i] = 0;
      X0[i] = 0;
      RM[i] = 0;
      for(j=0;j<3;++j) materialfraction[j][i] = 0;
      for(j=0;j<8;++j) layerfraction[j][i] = 0;
    }
  for(i=0;i<16;++i)
    for(j=0;j<8;++j)
      for(k=0;k<12;++k)
        OffSatu[i][j][k] = 0;
}

void FullShowerDevelopmentDescription::GetTrajectorySegment(double *pp, double *vv, double *ppstart, double *ppend)
{
  int i,j;
  double lambda;
  double pp0[3][3];
  double pp1[3][3];
  double segmentlength[3];
  // looking for the best axis
  int idir = -1;
  double segmentlengthmin = 99999999;

  for(i=0;i<3;++i)
    {
      if(vv[i]==0) continue;
      lambda = (m_fsgm->FSGM_CalSafeBoundaries[i][0]-pp[i])/vv[i];
      for(j=0;j<3;++j)
        {
          if(j!=i) pp0[i][j] = pp[j]+lambda*vv[j];
          else pp0[i][j] = m_fsgm->FSGM_CalSafeBoundaries[i][0];
        }
      lambda = (m_fsgm->FSGM_CalSafeBoundaries[i][1]-pp[i])/vv[i];
      for(j=0;j<3;++j)
        {
          if(j!=i) pp1[i][j] = pp[j]+lambda*vv[j];
          else pp1[i][j] = m_fsgm->FSGM_CalSafeBoundaries[i][1];
        }
      segmentlength[i] = 0;
      for(j=0;j<3;++j) segmentlength[i] += (pp0[i][j]-pp1[i][j])*(pp0[i][j]-pp1[i][j]);
      segmentlength[i] = sqrt(segmentlength[i]);
      if(segmentlength[i]<segmentlengthmin)
        {
          segmentlengthmin = segmentlength[i];
          idir = i;
        }
    }

  double checkdirection = (pp1[idir][0]-pp0[idir][0])*vv[0]+(pp1[idir][1]-pp0[idir][1])*vv[1]+(pp1[idir][2]-pp0[idir][2])*vv[2];
  
  if(checkdirection>0)
    {
      for(i=0;i<3;++i)
        {
          ppstart[i] = pp0[idir][i];
          ppend[i] = pp1[idir][i];
        }
    }
  else
    {
      for(i=0;i<3;++i)
        {
          ppstart[i] = pp1[idir][i];
          ppend[i] = pp0[idir][i];
        }
    }
}

bool FullShowerDevelopmentDescription::Compute(double *pp, double *vv, double startx0_input, double x0maxshower_input)
{
  if(Type!=0)
    {
      //  FullShowerDevelopmentDescription WRONG TYPE
      return false;
    }

  int i,j;
  int whereincal[4];
  // whereincal[0] = material
  // whereincal[1] = layer
  // whereincal[2] = tower
  // whereincal[3] = column

  NStep = 0;
  x0maxshower = x0maxshower_input;
  startx0 = startx0_input;
  lastx0 = 0;
  crackmaxfrac = 0;
  x0crackmaxfrac = 0;
  if(x0maxshower<2)
    {
      // FullShowerDevelopmentDescription x0maxshower is too small
      x0maxshower = 2;
    }

  double ppc[3];
  double ppcc[3];
  double pp0[3];
  double pp1[3];
  double vv0[3];
  double vv1[3];
  double vv2[3];
  double vnorm = sqrt(vv[0]*vv[0]+vv[1]*vv[1]+vv[2]*vv[2]);
  if(vnorm==0) return false;

  vv2[0] = vv[0]/vnorm;
  vv2[1] = vv[1]/vnorm;
  vv2[2] = vv[2]/vnorm;

  GetTrajectorySegment(pp,vv2,pp0,pp1);

  vv0[0] = vv2[1];
  vv0[1] = -vv2[0];
  vv0[2] = 0;
  vnorm = sqrt(vv0[0]*vv0[0]+vv0[1]*vv0[1]+vv0[2]*vv0[2]);
  if(vnorm>0)
    {
      vv0[0] = vv0[0]/vnorm;
      vv0[1] = vv0[1]/vnorm;
      vv0[2] = vv0[2]/vnorm;
    }
  else
    {
      vv0[0] = 1;
      vv0[1] = 0;
      vv0[2] = 0;
    }
  vv1[0] = vv2[1]*vv0[2]-vv2[2]*vv0[1];
  vv1[1] = vv2[2]*vv0[0]-vv2[0]*vv0[2];
  vv1[2] = vv2[0]*vv0[1]-vv2[1]*vv0[0];

  double lambda = sqrt( (pp1[0]-pp0[0])*(pp1[0]-pp0[0])+(pp1[1]-pp0[1])*(pp1[1]-pp0[1])+(pp1[2]-pp0[2])*(pp1[2]-pp0[2]) );

  int nhist = (int)floor(lambda/ZStep);
  if(nhist>=FSDD_NSTEPS_MAX)
    {
      // FullShowerDevelopmentDescription::Compute : nhist >= FSDD_NSTEPS_MAX
      return false;
    }
  NStep = nhist;

  double etotdep;

  double z2X0mat;

  X0[0] = startx0_input;

  double radius, relradius;
  double effradius, releffradius;
  double radialprofile;

  double x0position = X0[0];
  double relx0position;

  totx0cal = 0;
  totx0crack = 0;
  for(j=0;j<8;++j) totx0lay[j] = 0;
  for(j=0;j<8;++j) posx0lay[j] = 0;

  for(j=0;j<NStep;++j)
    {
      if(j>0) X0[j] = dX0[j-1] + X0[j-1];

      lambda = ZStep*((double)j+0.5);
      ppc[0] = pp0[0]+vv2[0]*lambda;
      ppc[1] = pp0[1]+vv2[1]*lambda;
      ppc[2] = pp0[2]+vv2[2]*lambda;
      for(i=0;i<3;++i) materialfraction[i][j] = 0;
      for(i=0;i<8;++i) layerfraction[i][j] = 0;
      etotdep = 0;
      relx0position = x0position/x0maxshower;
      releffradius = wideningfactor*m_fsgm->GetEffectiveRadius(relx0position,RadialContainedFraction);
      RM[j] = releffradius;
      effradius = 1.5*FSDD_MOLRAD*releffradius;
      for(i=0;i<m_fsgm->FSGM_NCircle;++i)
        {
          ppcc[0] = ppc[0]+effradius*(m_fsgm->FSGM_XCircle[i]*vv0[0]+m_fsgm->FSGM_YCircle[i]*vv1[0]);
          ppcc[1] = ppc[1]+effradius*(m_fsgm->FSGM_XCircle[i]*vv0[1]+m_fsgm->FSGM_YCircle[i]*vv1[1]);
          ppcc[2] = ppc[2]+effradius*(m_fsgm->FSGM_XCircle[i]*vv0[2]+m_fsgm->FSGM_YCircle[i]*vv1[2]);
          m_fsgm->WhereInCal(ppcc,whereincal);
          radius = m_fsgm->FSGM_RCircle[i]*effradius;
          relradius = radius/FSDD_MOLRAD;
          radialprofile = m_fsgm->GetRadialProfile(relx0position,relradius/wideningfactor);
          etotdep += radialprofile;
          materialfraction[whereincal[0]][j] += radialprofile;
          if(whereincal[0]==1)
            {
              layerfraction[whereincal[1]][j] += radialprofile;
            }
        }
      if(etotdep>0)
        {
          for(i=0;i<3;++i)
            materialfraction[i][j] /= etotdep;
          for(i=0;i<8;++i)
            layerfraction[i][j] /= etotdep;
        }
      
      z2X0mat = materialfraction[1][j]/FSDD_XCSI+materialfraction[2][j]/FSDD_XCRK;
      
      dX0[j] = ZStep*z2X0mat;
      x0position = X0[j];
      totx0cal += materialfraction[1][j]*dX0[j];
      totx0crack += materialfraction[2][j]*dX0[j];
      for(i=0;i<8;++i) totx0lay[i] += layerfraction[i][j]*dX0[j];

      if(materialfraction[2][j]>crackmaxfrac)
        {
          x0crackmaxfrac = X0[j];
          crackmaxfrac = materialfraction[2][j];
        }
    }

  dX0[NStep] = 0;
  X0[NStep] = dX0[NStep-1] + X0[NStep-1];
  lastx0 = X0[NStep];
  materialfraction[0][NStep] = 1;
  for(i=1;i<3;++i) materialfraction[i][NStep] = 0;
  for(i=0;i<8;++i) layerfraction[i][NStep] = 0;

  RemoveEmptySteps();

  return true;
}

void FullShowerDevelopmentDescription::RemoveEmptySteps()
{
  int i,j;
  int iwo0 = 0;
  for(i=0;i<NStep;++i)
    {
      if(dX0[i]==0) continue;
      if(iwo0==i) continue;
      X0[iwo0] = X0[i];
      dX0[iwo0] = dX0[i];
      RM[iwo0] = RM[i];
      for(j=0;j<3;++j) materialfraction[j][iwo0] = materialfraction[j][i];
      for(j=0;j<8;++j) layerfraction[j][iwo0] = layerfraction[j][i];
      ++iwo0;
    }
  dX0[iwo0] = 0;
  X0[iwo0] = X0[NStep];
  for(j=0;j<3;++j) materialfraction[j][iwo0] = materialfraction[j][NStep];
  for(j=0;j<8;++j) layerfraction[j][iwo0] = layerfraction[j][NStep];
  NStep = iwo0;
}

bool FullShowerDevelopmentDescription::ConvertToFixedX0(double x0step, FullShowerDevelopmentDescription *shmm)
{
  if(shmm==NULL) return false;
  if(shmm->NStep==0) return false;

  int i,j,ix0;
  double check[FSDD_NSTEPS_MAX];
  int icheck[FSDD_NSTEPS_MAX];

  NStep = 0;
  x0maxshower = shmm->x0maxshower;
  RadialContainedFraction = shmm->RadialContainedFraction ;
  startx0 = shmm->startx0;
  lastx0 = shmm->lastx0;
  crackmaxfrac = 0;
  x0crackmaxfrac = 0;
  totx0cal = shmm->totx0cal;
  totx0crack = 0;
  for(i=0;i<8;++i) totx0lay[i] = shmm->totx0lay[i];
  for(i=0;i<8;++i) posx0lay[i] = -1;

  ZStep = -1;

  for(i=0;i<FSDD_NSTEPS_MAX;++i)
    {
      check[i] = 0;
      icheck[i] = 0;
      X0[i] = x0step*(double)i;
      dX0[i] = x0step;
      for(j=0;j<3;++j) materialfraction[j][i] = 0;
      for(j=0;j<8;++j) layerfraction[j][i] = 0;  
      RM[i] = m_fsgm->GetEffectiveRadius(X0[i]/x0maxshower,RadialContainedFraction);
    }

  NStep = (int)floor(shmm->X0[shmm->NStep]/x0step)+1;
  if(NStep>=FSDD_NSTEPS_MAX)
    {
      //  FullShowerDevelopmentDescription::Convert : NStep >= FSDD_NSTEPS_MAX
      return false;
    }

  int start_ix0 = (int)floor(shmm->X0[0]/x0step);
  for(i=0;i<start_ix0;++i)
    materialfraction[0][i] = 1;

  int dimm = 0;
  int lastdimm = 0;

  double last_materialfraction[3];
  double last_layerfraction[8];

  last_materialfraction[0] = 1;
  for(j=1;j<3;++j) last_materialfraction[j] = 0;
  for(j=0;j<8;++j) last_layerfraction[j] = 0;

  int imm = 0;
  for(ix0=start_ix0;ix0<NStep;++ix0)
    {
      lastdimm = dimm;
      dimm = 0;
      while(shmm->X0[imm]<X0[ix0+1] && imm<shmm->NStep)
        {
          if(dimm==0)
            {
              check[ix0] += shmm->X0[imm]-X0[ix0];
              for(j=0;j<3;++j) materialfraction[j][ix0] += (shmm->X0[imm]-X0[ix0])*last_materialfraction[j];
              for(j=0;j<8;++j) layerfraction[j][ix0] += (shmm->X0[imm]-X0[ix0])*last_layerfraction[j];
            }
          icheck[ix0] += 1;
          if(shmm->X0[imm+1]<X0[ix0+1])
            {
              check[ix0] += shmm->dX0[imm];
              for(j=0;j<3;++j) materialfraction[j][ix0] += shmm->dX0[imm]*shmm->materialfraction[j][imm];
              for(j=0;j<8;++j) layerfraction[j][ix0] += shmm->dX0[imm]*shmm->layerfraction[j][imm];
            }
          else
            {
              check[ix0] += (X0[ix0+1]-shmm->X0[imm]);
              for(j=0;j<3;++j) materialfraction[j][ix0] += (X0[ix0+1]-shmm->X0[imm])*shmm->materialfraction[j][imm];
              for(j=0;j<8;++j) layerfraction[j][ix0] += (X0[ix0+1]-shmm->X0[imm])*shmm->layerfraction[j][imm];
            }
          for(j=0;j<3;++j) last_materialfraction[j] = shmm->materialfraction[j][imm];
          for(j=0;j<8;++j) last_layerfraction[j] = shmm->layerfraction[j][imm];
          ++imm;
          ++dimm;
        }
      if(dimm==0 && ix0<NStep-1)
        {
           for(j=0;j<3;++j) materialfraction[j][ix0] = last_materialfraction[j];
           for(j=0;j<8;++j) layerfraction[j][ix0] = last_layerfraction[j];
        }
    }

  ix0 = NStep-1;
  if(icheck[ix0]==0)
    {
      check[ix0] = shmm->X0[shmm->NStep]-X0[ix0];
      for(j=0;j<3;++j) materialfraction[j][ix0] = (shmm->X0[shmm->NStep]-X0[ix0])*last_materialfraction[j];
      for(j=0;j<8;++j) layerfraction[j][ix0] = (shmm->X0[shmm->NStep]-X0[ix0])*last_layerfraction[j];
    }
  check[ix0] += (X0[NStep]-shmm->X0[shmm->NStep]);
  for(j=0;j<3;++j) materialfraction[j][ix0] += (X0[NStep]-shmm->X0[shmm->NStep])*shmm->materialfraction[j][shmm->NStep];
  for(j=0;j<8;++j) layerfraction[j][ix0] += (X0[NStep]-shmm->X0[shmm->NStep])*shmm->layerfraction[j][shmm->NStep];
  
  for(i=NStep;i<FSDD_NSTEPS_MAX;++i)
    materialfraction[0][i] = 1;

  for(ix0=start_ix0;ix0<NStep;++ix0)
    {
      if(check[ix0]==0) continue;
      if(fabs(check[ix0]-x0step)>0.0001)
        {
          // FullShowerDevelopmentDescription::Convert PROBLEM DURING CONVERSION
          return false;
        }
      for(j=0;j<3;++j) materialfraction[j][ix0] /= x0step;
      for(j=0;j<8;++j) layerfraction[j][ix0] /= x0step;
    }

  double meanpos,meanweight;
  for(j=0;j<8;++j)
    {
      meanpos = 0;
      meanweight = 0;
      for(ix0=0;ix0<NStep;++ix0)
        {
          if(layerfraction[j][ix0]>0)
            {
              meanpos += layerfraction[j][ix0]*X0[ix0];
              meanweight += layerfraction[j][ix0];
            }
        }
      if(meanweight==0)
        posx0lay[j] = -1;
      else
        posx0lay[j] = meanpos/meanweight;
    }

  return true;
}

void FullShowerDevelopmentDescription::SetWideningFactor(double widfact)
{
  wideningfactor = widfact;
}

FullShowerDevelopmentDescriptionManager::FullShowerDevelopmentDescriptionManager(IGlastDetSvc *m_detSvc_input, int nxmax_input, double xmax0, double dxmax, double zstep_input, double radialcontainedfraction_input, double x0step)
  :m_detSvc(m_detSvc_input)
{

  // Initialize geometry
  m_fsgm = new FullShowerGeometryManager(m_detSvc);

  NDevelopment = 0;
  mintotx0cal = 0;
  maxtotx0cal = 0;
  int i;

  for(i=0;i<8;++i)
    {
      meantotx0lay[i] = 0;;
      meanposx0lay[i] = -1;
    }
  

  for(i=0;i<FSDD_NMAX;++i)
    {
      FSDDMM[i] = NULL;
      FSDDX0[i] = NULL;
    }
  CurrentFSDD = NULL;

  int nxmax = nxmax_input;
  if(nxmax>FSDD_NMAX-1)
    {
      // FullShowerDevelopmentDescriptionManager nxmax>FSDD_NMAX-1
      nxmax = FSDD_NMAX-1;
    }

  NDevelopment = nxmax;
  DXMax = dxmax;
  ZStep = zstep_input;;
  RadialContainedFraction = radialcontainedfraction_input;
  X0Step = x0step;
  
  for(i=0;i<=NDevelopment;++i)
    {
      XMax[i] = xmax0+DXMax*(double)i;
      FSDDMM[i] = new FullShowerDevelopmentDescription(m_fsgm,0,ZStep,RadialContainedFraction);
      FSDDX0[i] = new FullShowerDevelopmentDescription(m_fsgm,1,X0Step,RadialContainedFraction);
    }
  CurrentFSDD = new FullShowerDevelopmentDescription(m_fsgm,1,X0Step,RadialContainedFraction);
}

FullShowerDevelopmentDescriptionManager::~FullShowerDevelopmentDescriptionManager()
{
  int i;
  for(i=0;i<FSDD_NMAX;++i)
    {
      if(FSDDMM[i]!=NULL) delete FSDDMM[i];
      if(FSDDX0[i]!=NULL) delete FSDDX0[i];
    }
  if(CurrentFSDD!=NULL) delete CurrentFSDD;
}

bool FullShowerDevelopmentDescriptionManager::Compute(double *pp, double *vv, double startx0_input)
{
  int i,j,k,l;
  mintotx0cal = 99999999;
  maxtotx0cal = -99999999;

  for(l=0;l<=NDevelopment;++l)
    {
      for(i=0;i<16;++i)
        for(j=0;j<8;++j)
          for(k=0;k<12;++k)
            {
              FSDDMM[l]->OffSatu[i][j][k] = OffSatu[i][j][k];
              FSDDX0[l]->OffSatu[i][j][k] = OffSatu[i][j][k];
            }
    }
  
  for(i=0;i<=NDevelopment;++i)
    {
      if(!FSDDMM[i]->Compute(pp,vv,startx0_input,XMax[i])) return false;
      if(!FSDDX0[i]->ConvertToFixedX0(X0Step,FSDDMM[i])) return false;
      if(FSDDX0[i]->totx0cal<mintotx0cal) mintotx0cal = FSDDX0[i]->totx0cal;
      if(FSDDX0[i]->totx0cal>maxtotx0cal) maxtotx0cal = FSDDX0[i]->totx0cal;
    }
  
  double meanw = 0;
  for(j=0;j<8;++j)
    {
      meantotx0lay[j] = 0;;
      meanposx0lay[j] = 0;
      meanw = 0;
      for(i=0;i<=NDevelopment;++i)
        {
          meantotx0lay[j] += FSDDX0[i]->totx0lay[j];
          if(FSDDX0[i]->posx0lay[j]>-1)
            {
              meanw += 1.;
              meanposx0lay[j] += FSDDX0[i]->posx0lay[j];
            }
        }
      meantotx0lay[j] /= (double)(NDevelopment+1);
      if(meanw>0)
        meanposx0lay[j] /= meanw;
      else
        meanposx0lay[j] = -1;
    }
  
  return true;
}

void FullShowerDevelopmentDescriptionManager::FillCurrentFSDD(double showerxmax)
{  
  int ish;
  double interpol;

  if(showerxmax<=XMax[0])
    {
      ish = 0;
      interpol = 0;
    }
  else if(showerxmax>=XMax[NDevelopment])
    {
      ish = NDevelopment-1;
      interpol = 1;
    }
  else
    {
      ish = (int)floor((showerxmax-XMax[0])/DXMax);
      interpol = (showerxmax-XMax[ish])/DXMax;
    }

  int i,j;
  CurrentFSDD->NStep = FSDDX0[ish]->NStep;
  if(FSDDX0[ish+1]->NStep>CurrentFSDD->NStep)
    CurrentFSDD->NStep = FSDDX0[ish+1]->NStep;
  CurrentFSDD->x0maxshower = showerxmax;
  CurrentFSDD->crackmaxfrac = (1-interpol)*FSDDX0[ish]->crackmaxfrac+interpol*FSDDX0[ish+1]->crackmaxfrac;
  CurrentFSDD->x0crackmaxfrac = (1-interpol)*FSDDX0[ish]->x0crackmaxfrac+interpol*FSDDX0[ish+1]->x0crackmaxfrac;
  CurrentFSDD->startx0 = (1-interpol)*FSDDX0[ish]->startx0+interpol*FSDDX0[ish+1]->startx0;
  CurrentFSDD->lastx0 = (1-interpol)*FSDDX0[ish]->lastx0+interpol*FSDDX0[ish+1]->lastx0;
  CurrentFSDD->totx0cal = (1-interpol)*FSDDX0[ish]->totx0cal+interpol*FSDDX0[ish+1]->totx0cal;
  CurrentFSDD->totx0crack = (1-interpol)*FSDDX0[ish]->totx0crack+interpol*FSDDX0[ish+1]->totx0crack;
  for(i=0;i<8;++i) CurrentFSDD->totx0lay[i] = (1-interpol)*FSDDX0[ish]->totx0lay[i]+interpol*FSDDX0[ish+1]->totx0lay[i];
  
  for(i=0;i<=CurrentFSDD->NStep;++i)
    {
      CurrentFSDD->dX0[i] = FSDDX0[ish]->dX0[i];
      CurrentFSDD->X0[i] = FSDDX0[ish]->X0[i];
      CurrentFSDD->RM[i] = FSDDX0[ish]->RM[i];
      for(j=0;j<3;++j)
        CurrentFSDD->materialfraction[j][i] = (1-interpol)*FSDDX0[ish]->materialfraction[j][i]
          +interpol*FSDDX0[ish+1]->materialfraction[j][i];
      for(j=0;j<8;++j)
        CurrentFSDD->layerfraction[j][i] = (1-interpol)*FSDDX0[ish]->layerfraction[j][i]
          +interpol*FSDDX0[ish+1]->layerfraction[j][i];
    }
/*   for(i=CurrentFSDD->NStep;i<FSDD_NSTEPS_MAX;++i) */
/*     { */
/*       CurrentFSDD->dX0[i] = FSDDX0[ish]->dX0[i]; */
/*       CurrentFSDD->X0[i] = FSDDX0[ish]->X0[i]; */
/*       CurrentFSDD->RM[i] = FSDDX0[ish]->RM[i]; */
/*       CurrentFSDD->materialfraction[0][i] = 1; */
/*       for(j=1;j<3;++j) CurrentFSDD->materialfraction[j][i] = 0; */
/*       for(j=0;j<8;++j) CurrentFSDD->layerfraction[j][i] = 0; */
/*     } */
  return;
}

void FullShowerDevelopmentDescriptionManager::SetWideningFactor(double widfact)
{  

  int i;
  for(i=0;i<=NDevelopment;++i)
    {
      FSDDMM[i]->SetWideningFactor(widfact);
      FSDDX0[i]->SetWideningFactor(widfact);
    }
  CurrentFSDD->SetWideningFactor(widfact);
}

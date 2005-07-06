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

#include "TMath.h"

/**   
* @class FullShowerDevelopmentDescription
* @author Philippe Bruel
*
* Tool that describes the shower developement in the calorimeter given
* the length in X0 seen in the tracker and the position of the shower maximum
*
* $Header$
*/

FullShowerDevelopmentDescription::FullShowerDevelopmentDescription(IGlastDetSvc *m_detSvc_input, int type_input, double zstep_input, double radialcontainedfraction_input)
  :m_detSvc(m_detSvc_input)
{
  
  Initialize();
  Reset();
  FSDD_XCircle = NULL;
  FSDD_YCircle = NULL;
  FSDD_RCircle = NULL;
  Type = type_input;
  if(Type==0) PrepXYPoints();
  ZStep = zstep_input;
  RadialContainedFraction = radialcontainedfraction_input;
}

FullShowerDevelopmentDescription::~FullShowerDevelopmentDescription()
{
  
}

void FullShowerDevelopmentDescription::Initialize()
{

  // Geometry

  FSDD_CsIHeight = 19.9;
  m_detSvc->getNumericConstByName(std::string("CsIHeight"),&FSDD_CsIHeight);

  FSDD_CsIWidth = 26.7;
  m_detSvc->getNumericConstByName(std::string("CsIWidth"),&FSDD_CsIWidth);

  FSDD_CsILength = 326.0;
  m_detSvc->getNumericConstByName(std::string("CsILength"),&FSDD_CsILength);

  FSDD_cellVertPitch = 21.35;
  m_detSvc->getNumericConstByName(std::string("cellVertPitch"),&FSDD_cellVertPitch);

  FSDD_cellHorPitch = 27.84;
  m_detSvc->getNumericConstByName(std::string("cellHorPitch"),&FSDD_cellHorPitch);

  FSDD_towerPitch = 374.5;
  m_detSvc->getNumericConstByName(std::string("towerPitch"),&FSDD_towerPitch);

  // FSDD_crackHalfWidth is the average of X and Y layers crackHalfWidths
  FSDD_crackHalfWidth = ((FSDD_towerPitch-FSDD_CsILength)/2 + (FSDD_towerPitch-12*FSDD_cellHorPitch)/2)/2;;

  // Retrieving the top position of the calorimeter as in TkrUtil/*/src/TkrGeometrySvc.cxx
  FSDD_calZTop = -47.395;

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
  HepTransform3D transfTop;
  int count;
  for (count=0;count<3;++count) {
    topLayerId.append(0);
    if((sc = m_detSvc->getTransform3DByID(topLayerId,&transfTop)).isSuccess()) break;
  }
  if(sc.isSuccess()) 
    {
      FSDD_calZTop = (transfTop.getTranslation()).z()+ FSDD_cellVertPitch/2;
    }

  FSDD_calZBot = FSDD_calZTop-8*FSDD_cellVertPitch;

  double safedistance = 50;
  FSDD_CalSafeBoundaries[0][0] = -2*FSDD_towerPitch-safedistance;
  FSDD_CalSafeBoundaries[0][1] = +2*FSDD_towerPitch+safedistance;
  FSDD_CalSafeBoundaries[1][0] = -2*FSDD_towerPitch-safedistance;
  FSDD_CalSafeBoundaries[1][1] = +2*FSDD_towerPitch+safedistance;
  FSDD_CalSafeBoundaries[2][0] = FSDD_calZBot-safedistance;
  FSDD_CalSafeBoundaries[2][1] = FSDD_calZTop+safedistance;

  // Radiation lengths
  FSDD_MOLRAD = 35.;
  FSDD_gCSI = 4.51;
  FSDD_gCRK = 0.607287449392712508; // 10.*2.7/(2.*22.23) // 10mm
  FSDD_XCSI = 18.5*(FSDD_cellVertPitch/FSDD_CsIHeight);
  FSDD_XCRK = 395.364666666666722; // 24.01/FSDD_gCRK*10. // *10 to get it in mm

  // Parameters used for the radial shower development description
  FSDD_CORE0 = 0.0211192;
  FSDD_CORE1 = 0.114962;
  FSDD_TAIL0 = 0.450839;
  FSDD_TAIL1 = 0.840738;
  FSDD_TAIL2 = -1.87526;
  FSDD_TAIL3 = 0.643606;
  FSDD_PCT0 = 2.51719;
  FSDD_PCT1 = 0.554606;
  FSDD_PCT2 = 0.851856;

}

void FullShowerDevelopmentDescription::Reset()
{

  int i,j;
  NStep = 0;
  ZStep = 0;
  x0maxshower = 0;
  crackmaxfrac = 0;
  x0crackmaxfrac = 0;
  firstx0 = 0;
  lastx0 = 0;
  totx0cal = 0;
  totx0crack = 0;
  for(i=0;i<8;++i) x0lay[i] = 0;
  for(i=0;i<FSDD_NSTEPS_MAX;++i)
    {
      dX0[i] = 0;
      X0[i] = 0;
      RM[i] = 0;
      for(j=0;j<3;++j) materialfraction[j][i] = 0;
      for(j=0;j<8;++j) layerfraction[j][i] = 0;
    }
}

void FullShowerDevelopmentDescription::PrepXYPoints()
{
  if(FSDD_XCircle==NULL) FSDD_XCircle = new double[FSDD_NPOINTS_MAX];
  if(FSDD_YCircle==NULL) FSDD_YCircle = new double[FSDD_NPOINTS_MAX];
  if(FSDD_RCircle==NULL) FSDD_RCircle = new double[FSDD_NPOINTS_MAX];

  int n = (int)sqrt(1/3.2*(double)FSDD_NPOINTS_MAX);
  int i,j;
  double step = 1/(double)(n);
  double x,y;
  FSDD_NCircle = 0;
  for(i=0;i<=2*n;++i)
    {
      x = -1+step*(double)i;
      for(j=0;j<=2*n;++j)
	{
	  y = -1+step*(double)j;
	  if(sqrt(x*x+y*y)>1) continue;
	  FSDD_XCircle[FSDD_NCircle] = x;
	  FSDD_YCircle[FSDD_NCircle] = y;
	  FSDD_RCircle[FSDD_NCircle] = sqrt(x*x+y*y);
	  ++FSDD_NCircle;
	}
    }
}

double FullShowerDevelopmentDescription::RCore(double x)
{
  double xx = x;
  if(xx<.25) xx = .25;
  if(xx>2.5) xx = 2.5;
  double val = FSDD_CORE0+FSDD_CORE1*xx;
  return val;
}

double FullShowerDevelopmentDescription::RTail(double x)
{
  double xx = x;
  if(xx<.25) xx = .25;
  if(xx>2.5) xx = 2.5;
  return FSDD_TAIL0*(exp(FSDD_TAIL2*(xx-FSDD_TAIL1))+exp(FSDD_TAIL3*(xx-FSDD_TAIL1)));
}

double FullShowerDevelopmentDescription::PCore(double x)
{
  double xx = x;
  if(xx<.25) xx = .25;
  if(xx>2.5) xx = 2.5;
  double arg = FSDD_PCT0*exp((FSDD_PCT1-xx)/FSDD_PCT2-exp((FSDD_PCT1-xx)/FSDD_PCT2));
  if(arg<0) arg = 0;
  if(arg>1) arg = 1;
  return arg;
}

double FullShowerDevelopmentDescription::RadialProfile(double r, double x)
{
  double rcore2 = RCore(x);
  rcore2 = rcore2*rcore2;
  double rtail2 = RTail(x);
  rtail2 = rtail2*rtail2;
  double pcore = PCore(x);
  return pcore*(2*rcore2/(r*r+rcore2)/(r*r+rcore2))
    + (1-pcore)*(2*rtail2/(r*r+rtail2)/(r*r+rtail2));
}

double FullShowerDevelopmentDescription::EffectiveRadius(double z, double radialcontainedfraction)
{
  double rcore2 = RCore(z);
  rcore2 = rcore2*rcore2;
  double rtail2 = RTail(z);
  rtail2 = rtail2*rtail2;
  double pcore = PCore(z);
  
  double a = (1-radialcontainedfraction)/rcore2/rtail2;
  double b = (1-radialcontainedfraction)*(1/rcore2+1/rtail2) - pcore/rtail2 - (1-pcore)/rcore2;
  double c = -radialcontainedfraction;
  double delta = b*b-4*a*c;
  double r1 = (-b+sqrt(delta))/2/a;
  return sqrt(r1);
}

void FullShowerDevelopmentDescription::WhereInCal(double *xyz, int *whereincal)
{
  whereincal[0] = 0;
  whereincal[1] = -1;
  if(xyz[2]>FSDD_calZTop) return;
  if(xyz[2]<FSDD_calZBot) return;
  if(fabs(xyz[0])>2*FSDD_towerPitch) return;
  if(fabs(xyz[1])>2*FSDD_towerPitch) return;
  double xmod = xyz[0]+2*FSDD_towerPitch-FSDD_towerPitch*floor((xyz[0]+2*FSDD_towerPitch)/FSDD_towerPitch);
  if(FSDD_towerPitch-xmod<xmod) xmod = FSDD_towerPitch-xmod;
  double ymod = xyz[1]+2*FSDD_towerPitch-FSDD_towerPitch*floor((xyz[1]+2*FSDD_towerPitch)/FSDD_towerPitch);
  if(FSDD_towerPitch-ymod<ymod) ymod = FSDD_towerPitch-ymod;
  if(xmod<FSDD_crackHalfWidth || ymod<FSDD_crackHalfWidth)
    {
      whereincal[0] = 2;
      return;
    }
  int ilayer = (int)floor( (FSDD_calZTop-xyz[2])/FSDD_cellVertPitch );
  double xy = ymod;
  if(ilayer%2==1) xy = xmod;
  xy -= (FSDD_towerPitch/2-6*FSDD_cellHorPitch);
  double xymod = xy-FSDD_cellHorPitch*floor(xy/FSDD_cellHorPitch);
  if(FSDD_cellHorPitch-xymod<xymod) xymod = FSDD_cellHorPitch-xymod;
  if(xymod<(FSDD_cellHorPitch-FSDD_CsIWidth)/2) return;
/*   double zlayer = (FSDD_calZTop-xyz[2])-FSDD_cellVertPitch*(double)ilayer; */
/*   if(FSDD_cellVertPitch-zlayer>zlayer) zlayer = FSDD_cellVertPitch-zlayer; */
/*   if(zlayer<(FSDD_cellVertPitch-FSDD_CsIHeight)/2) return; */
  whereincal[0] = 1;
  whereincal[1] = ilayer;
}

void FullShowerDevelopmentDescription::Print(int i)
{
  printf("FullShowerDevelopmentDescription::Print  %d (RM = %4.2f) : %3.2f %3.2f %3.2f | %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f (dX0 = %f X0 = %f) \n",
	 i,RM[i],materialfraction[0][i],materialfraction[1][i],materialfraction[2][i],
	 layerfraction[0][i],layerfraction[1][i],layerfraction[2][i],layerfraction[3][i],
	 layerfraction[4][i],layerfraction[5][i],layerfraction[6][i],layerfraction[7][i],
	 dX0[i],X0[i]);  
}

void FullShowerDevelopmentDescription::PrintAll()
{
  int i;
  for(i=0;i<=NStep;++i) Print(i);
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
      lambda = (FSDD_CalSafeBoundaries[i][0]-pp[i])/vv[i];
      for(j=0;j<3;++j)
	{
	  if(j!=i) pp0[i][j] = pp[j]+lambda*vv[j];
	  else pp0[i][j] = FSDD_CalSafeBoundaries[i][0];
	}
      lambda = (FSDD_CalSafeBoundaries[i][1]-pp[i])/vv[i];
      for(j=0;j<3;++j)
	{
	  if(j!=i) pp1[i][j] = pp[j]+lambda*vv[j];
	  else pp1[i][j] = FSDD_CalSafeBoundaries[i][1];
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

bool FullShowerDevelopmentDescription::Compute(double *pp, double *vv, double startx0, double x0maxshower_input)
{
  if(Type!=0)
    {
      //      printf("FullShowerDevelopmentDescription WRONG TYPE !!!!!\n");
      return false;
    }

  int i,j;
  int whereincal[2];

  NStep = 0;
  x0maxshower = x0maxshower_input;
  firstx0 = startx0;
  lastx0 = 0;
  crackmaxfrac = 0;
  x0crackmaxfrac = 0;
  if(x0maxshower<2)
    {
      //      printf("FullShowerDevelopmentDescription x0maxshower is too small : %f !!! ... forcing it to be = 2\n",x0maxshower);
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
      //      printf("FullShowerDevelopmentDescription::Compute : nhist >= FSDD_NSTEPS_MAX : %d !!!!!!!!!!!!!!!!!!!!!!\n",nhist);
      return false;
    }
  NStep = nhist;

  double etotdep;

  double z2X0mat;

  X0[0] = startx0;

  double radius, relradius;
  double effradius, releffradius;
  double radialprofile;

  double x0position = X0[0];
  double relx0position;

  totx0cal = 0;
  totx0crack = 0;
  for(j=0;j<8;++j) x0lay[j] = 0;

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
      releffradius = EffectiveRadius(relx0position,RadialContainedFraction);
      RM[j] = releffradius;
      effradius = 1.5*FSDD_MOLRAD*releffradius;
      for(i=0;i<FSDD_NCircle;++i)
	{
	  ppcc[0] = ppc[0]+effradius*FSDD_XCircle[i]*vv0[0]+effradius*FSDD_YCircle[i]*vv1[0];
	  ppcc[1] = ppc[1]+effradius*FSDD_XCircle[i]*vv0[1]+effradius*FSDD_YCircle[i]*vv1[1];
	  ppcc[2] = ppc[2]+effradius*FSDD_XCircle[i]*vv0[2]+effradius*FSDD_YCircle[i]*vv1[2];
	  WhereInCal(ppcc,whereincal);
	  radius = FSDD_RCircle[i]*effradius;
	  relradius = radius/35.;
	  radialprofile = RadialProfile(relradius,relx0position);
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
	    {
	      layerfraction[i][j] /= etotdep;
	    }
	}
      
      z2X0mat = materialfraction[1][j]/FSDD_XCSI+materialfraction[2][j]/FSDD_XCRK;
      
      dX0[j] = ZStep*z2X0mat;
      x0position = X0[j];
      totx0cal += materialfraction[1][j]*dX0[j];
      totx0crack += materialfraction[2][j]*dX0[j];
      for(i=0;i<8;++i) x0lay[i] += layerfraction[i][j]*dX0[j];

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
  firstx0 = shmm->firstx0;
  lastx0 = shmm->lastx0;
  crackmaxfrac = 0;
  x0crackmaxfrac = 0;
  totx0cal = shmm->totx0cal;
  totx0crack = 0;
  for(i=0;i<8;++i) x0lay[i] = shmm->x0lay[i];

  ZStep = -1;

  for(i=0;i<FSDD_NSTEPS_MAX;++i)
    {
      check[i] = 0;
      icheck[i] = 0;
      X0[i] = x0step*(double)i;
      dX0[i] = x0step;
      for(j=0;j<3;++j) materialfraction[j][i] = 0;
      for(j=0;j<8;++j) layerfraction[j][i] = 0;  
      RM[i] = EffectiveRadius(X0[i]/x0maxshower,RadialContainedFraction);
    }

  NStep = (int)floor(shmm->X0[shmm->NStep]/x0step)+1;
  if(NStep>=FSDD_NSTEPS_MAX)
    {
      //      printf("FullShowerDevelopmentDescription::Convert !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! : NStep >= FSDD_NSTEPS_MAX\n");
      return false;
    }

  int start_ix0 = (int)floor(shmm->X0[0]/x0step);
  for(i=0;i<start_ix0;++i)
    materialfraction[0][i] = 1;

  double coherencetest;
  int dimm = 0;
  int lastdimm = 0;
  double firstx0;

  double last_materialfraction[3];
  double last_layerfraction[8];
  double last_x0;

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
	  //	  printf("FullShowerDevelopmentDescription::Convert PROBLEM DURING CONVERSION : check[%d] = %f !!!!!!!!!!!!!!!!\n",ix0,check[ix0]);
	  return false;
	}
      for(j=0;j<3;++j) materialfraction[j][ix0] /= x0step;
      for(j=0;j<8;++j) layerfraction[j][ix0] /= x0step;
    }

  return true;
}

FullShowerDevelopmentDescriptionManager::FullShowerDevelopmentDescriptionManager(IGlastDetSvc *m_detSvc_input, int nxmax_input, double xmax0, double dxmax, double zstep_input, double radialcontainedfraction_input, double x0step)
  :m_detSvc(m_detSvc_input)
{
  NDevelopment = 0;
  maxtotx0cal = 0;
  int i;
  for(i=0;i<FSDD_NMAX;++i)
    {
      FSDDMM[i] = NULL;
      FSDDX0[i] = NULL;
    }
  CurrentFSDD = NULL;

  int nmax = nmax_input;
  if(nxmax>FSDD_NMAX-1)
    {
      //      printf("FullShowerDevelopmentDescriptionManager nxmax>FSDD_NMAX-1 !!!!!!!!!!!!!!!!\n");
      nmax = FSDD_NMAX-1;
    }

  NDevelopment = nxmax;
  DXMax = dxmax;
  ZStep = zstep_input;;
  RadialContainedFraction = radialcontainedfraction_input;
  X0Step = x0step;
  
  for(i=0;i<=NDevelopment;++i)
    {
      XMax[i] = xmax0+DXMax*(double)i;
      FSDDMM[i] = new FullShowerDevelopmentDescription(m_detSvc,0,ZStep,RadialContainedFraction);
      FSDDX0[i] = new FullShowerDevelopmentDescription(m_detSvc,1,X0Step,0);
    }
  CurrentFSDD = new FullShowerDevelopmentDescription(m_detSvc,1,X0Step,0);
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

bool FullShowerDevelopmentDescriptionManager::Compute(double *pp, double *vv, double startx0)
{
  int i;
  maxtotx0cal = 0;
  for(i=0;i<=NDevelopment;++i)
    {
      if(!FSDDMM[i]->Compute(pp,vv,startx0,XMax[i])) return false;
      if(!FSDDX0[i]->ConvertToFixedX0(X0Step,FSDDMM[i])) return false;
      if(FSDDX0[i]->totx0cal>maxtotx0cal) maxtotx0cal = FSDDX0[i]->totx0cal;
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
  CurrentFSDD->firstx0 = (1-interpol)*FSDDX0[ish]->firstx0+interpol*FSDDX0[ish+1]->firstx0;
  CurrentFSDD->lastx0 = (1-interpol)*FSDDX0[ish]->lastx0+interpol*FSDDX0[ish+1]->lastx0;
  CurrentFSDD->totx0cal = (1-interpol)*FSDDX0[ish]->totx0cal+interpol*FSDDX0[ish+1]->totx0cal;
  CurrentFSDD->totx0crack = (1-interpol)*FSDDX0[ish]->totx0crack+interpol*FSDDX0[ish+1]->totx0crack;
  for(i=0;i<8;++i) CurrentFSDD->x0lay[i] = (1-interpol)*FSDDX0[ish]->x0lay[i]+interpol*FSDDX0[ish+1]->x0lay[i];
  
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

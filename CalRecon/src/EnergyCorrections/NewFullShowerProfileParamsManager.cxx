/**   
*
* THIS CODE IS STILL NOT FULLY DOCUMENTED !!!!!!!!!!!!!!!!!!!!
* DESCRIPTION COMMENTS WILL BE ADDED SOON.
* Philippe Bruel
*
*/

#include "NewFullShowerProfileParamsManager.h"

#include "TMath.h"

/**   
* @class NewFullShowerProfileParamsManager
* @author Philippe Bruel
*
* Class handling the energy parameterization of the shower profile parameters alpha and tmax;
*
* $Header$
*/

double mypol21(double x, Double_t *par)
{
  if(x>par[2]) return par[0]+par[1]*x;
  double x0 = par[2];
  double g1 = par[3];
  double b1 = par[1]-2*g1*x0;
  double a1 = par[0]+par[1]*x0-b1*x0-g1*x0*x0;
  return a1+b1*x+g1*x*x;
}

NewFullShowerProfileParamsManager::NewFullShowerProfileParamsManager(int optmodprofin)
{
  optmodprof = optmodprofin;
  LNA0 = 2.532670;
  LNA1 = 0.615655;
  SLNA0 = 2.946404;
  SLNA1 = 1.017415;
  LNT0 = 3.329578;
  SLNT0 = 4.043513;
  SLNT1 = 1.577324;
  RHO0 = 0.489448;
  RHO1 = -0.076717;
  RHO2 = 14.326741;
  RELERR0 = 0.211965;
  RELERR1 = -0.393502;
  A0 = 2.162854;
  A1 = 0.592538;
  T0 = 3.112031;
  T1 = 1.090746;

  THETACOR = atan(-0.5);
  CTHETACOR = cos(THETACOR);
  STHETACOR = sin(THETACOR);
  
  XCORMEAN = 0;
  YCORMEAN = 0;
  XCORSIGMA = 0;
  YCORSIGMA = 0;

  if(optmodprof==0)
    {
      // no modprof FORPAPER nomodprof newxycor with nbinx0=10
      XCORMPAR[0] = 1.058104e+00;
      XCORMPAR[1] = 4.068198e-01;
      XCORMPAR[2] = -3.613821e-02;
      XCORSPAR[0] = 1.941007e-01;
      XCORSPAR[1] = -2.882820e-02;
      XCORSPAR[2] = 2.138814e+00;
      XCORSPAR[3] = 2.454887e-02;
      YCORMPAR[0] = 1.760027e-02;
      YCORMPAR[1] = -1.695235e-01;
      YCORMPAR[2] = 1.471772e-02;
      YCORSPAR[0] = 5.013077e-02;
      YCORSPAR[1] = -7.635890e-03;
      YCORSPAR[2] = 1.961581e+00;
      YCORSPAR[3] = 8.104348e-03;
    }
  else
    {
      // with modprof
      XCORMPAR[0] = 1.084104e+00;
      XCORMPAR[1] = 3.906287e-01;
      XCORMPAR[2] = -3.330708e-02;
      XCORSPAR[0] = 1.954132e-01;
      XCORSPAR[1] = -2.950297e-02;
      XCORSPAR[2] = 1.973298e+00;
      XCORSPAR[3] = 2.940200e-02;
      YCORMPAR[0] = 2.660861e-02;
      YCORMPAR[1] = -1.753368e-01;
      YCORMPAR[2] = 1.579789e-02;
      YCORSPAR[0] = 5.158166e-02;
      YCORSPAR[1] = -8.028777e-03;
      YCORSPAR[2] = 1.837070e+00;
      YCORSPAR[3] = 1.118694e-02;
    }

  // superseeded by PAPER results from allshowh2errr(10) no modprof and NEW PARAM [0]+exp([1]+[2]*x)
  XYERRPAR[0] = 0.0;
  XYERRPAR[1] = -1.753371e+00;
  XYERRPAR[2] = -7.248170e-01;
  XYERRPAR[3] = 0;
  //
  XYERRPARL[0] = 2.1222;
  XYERRPARL[1] = -0.248913;
  //
  XYERRPARR[0] =  7.47350e-02;
  XYERRPARR[1] =  2.95238e-01;
  XYERRPARR[2] = -1.93662e+00;

  MODPROF00 = -0.293944;
  MODPROF01 = 0.201803;
  MODPROF02 = -0.0250792;
  
  MODPROF20 = 5.58211e-02;
  MODPROF21 =-2.16670e-02;
  MODPROF22 = 2.81842e-03;
  MODPROF23 = 7.13075e-01;
  MODPROF24 =-6.84458e-02;
  
  Reset();
}

NewFullShowerProfileParamsManager::~NewFullShowerProfileParamsManager()
{
  
}

void NewFullShowerProfileParamsManager::Reset()
{
  logenergy = -1;
  energy = -1;
  alpha = -1;
  lnalpha = -1;
  slnalpha = -1;
  tmax = -1;
  lntmax = -1;
  slntmax = -1;
  rho = -1;
  relerr = -1;
  relerrparl = -1;
  relerrparr = -1;
  beta = -1;
  logenergy = -1;
  xcor = -1;
  ycor = -1;

  XCORMEAN = -1;
  YCORMEAN = -1;
  XCORSIGMA = -1;
  YCORSIGMA = -1;

  modparam0 = -1;
  modparam1 = -1;
  modparam2 = -1;
}

bool NewFullShowerProfileParamsManager::Fill(double *par)
{
  logenergy = par[2];
  if(fabs(par[2])>307) {Reset(); return false;}

  energy = exp(log(10.)*logenergy);

  double loge = logenergy;
  if(loge>5) loge = 5;
  XCORMEAN = XCORMPAR[0]+XCORMPAR[1]*loge+XCORMPAR[2]*loge*loge;
  YCORMEAN = YCORMPAR[0]+YCORMPAR[1]*loge+YCORMPAR[2]*loge*loge;
  XCORSIGMA = mypol21(loge,XCORSPAR);
  YCORSIGMA = mypol21(loge,YCORSPAR);
  relerr = XYERRPAR[0]+exp(XYERRPAR[1]+loge*XYERRPAR[2]);
  relerrparl = XYERRPARL[0]+loge*XYERRPARL[1];
  relerrparr = XYERRPARR[0]+XYERRPARR[1]*exp(loge*XYERRPARR[2]);

  xcor = XCORMEAN+XCORSIGMA*par[0];
  ycor = YCORMEAN+YCORSIGMA*par[1];
  //
  double myloge = loge;
  if(myloge<0) myloge = 0;
  if(myloge>3.8) myloge = 3.8;
  modparam0 = MODPROF00+MODPROF01*myloge+MODPROF02*myloge*myloge;
  modparam1 = 0.9;
  if(myloge>MODPROF23)
    modparam2 = MODPROF20+MODPROF21*myloge+MODPROF22*myloge*myloge;
  else
    modparam2 = MODPROF20+MODPROF21*MODPROF23+MODPROF22*MODPROF23*MODPROF23+MODPROF24*(myloge-MODPROF23);
  //
  lnalpha = xcor*CTHETACOR + ycor*STHETACOR;
  if(lnalpha>709) {Reset(); return false;}

  alpha =  exp(lnalpha);
  beta = -xcor*STHETACOR + ycor*CTHETACOR;

  if(alpha<=1.) {Reset(); return false;}
  if(beta<=0.) {Reset(); return false;}

  tmax = (alpha-1)/beta;

  return true;
}

double NewFullShowerProfileParamsManager::GetChi2Contribution(double *par)
{
  return par[0]*par[0]+par[1]*par[1];
}

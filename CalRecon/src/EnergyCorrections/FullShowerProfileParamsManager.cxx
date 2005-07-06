/**   
*
* THIS CODE IS STILL NOT FULLY DOCUMENTED !!!!!!!!!!!!!!!!!!!!
* DESCRIPTION COMMENTS WILL BE ADDED SOON.
* Philippe Bruel
*
*/

#include "FullShowerProfileParamsManager.h"

#include "TMath.h"

/**   
* @class FullShowerProfileParamsManager
* @author Philippe Bruel
*
* Class handling the energy parameterization of the shower profile parameters alpha and tmax;
*
* $Header$
*/

FullShowerProfileParamsManager::FullShowerProfileParamsManager()
{
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

  Reset();
}

FullShowerProfileParamsManager::~FullShowerProfileParamsManager()
{
  
}

void FullShowerProfileParamsManager::Reset()
{
  energy = -1;
  alpha = -1;
  lnalpha = -1;
  slnalpha = -1;
  tmax = -1;
  lntmax = -1;
  slntmax = -1;
  rho = -1;
  relerr = -1;
  beta = -1;
}

bool FullShowerProfileParamsManager::Fill(double energy_input)
{
  if(energy_input<=0)
    {
      Reset();
      return false;
    }
  energy = energy_input;
  alpha = A1*log(energy)+A0;
  tmax = T1*log(energy)+T0;
  lnalpha = log(LNA1*log(energy)+LNA0);
  slnalpha = 1./(SLNA1*log(energy)+SLNA0);
  lntmax = log(log(energy)+LNT0);
  slntmax = 1./(SLNT1*log(energy)+SLNT0);
  rho = RHO0-RHO1*exp(-energy/RHO2);
  relerr = RELERR0*exp(RELERR1*log(energy));
  if(relerr>0.15) relerr = 0.15;
  beta = (alpha-1.)/tmax;
  return true;
}

double FullShowerProfileParamsManager::GetChi2Contribution(double *par)
{
  // par[0] = alpha_fit
  // par[1] = tmax_fit
  // par[2] = energy_fit
  if(!Fill(par[2])) return -1;
  if(par[0]<=1 || par[1]<=0) return -1;
  double argalpha = (log(par[0])-lnalpha)/slnalpha;
  double argtmax = (log(par[1])-lntmax)/slntmax;
  return (argalpha*argalpha + argtmax*argtmax -2*rho*argalpha*argtmax)/(1-rho*rho);
}

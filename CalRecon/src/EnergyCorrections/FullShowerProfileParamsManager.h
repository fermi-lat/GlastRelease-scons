/**   
*
* THIS CODE IS STILL NOT FULLY DOCUMENTED !!!!!!!!!!!!!!!!!!!!
* DESCRIPTION COMMENTS WILL BE ADDED SOON.
* Philippe Bruel
*
*/

/**   
* @class FullShowerProfileParamsManager
* @author Philippe Bruel
*
* Class handling the energy parameterization of the shower profile parameters alpha and tmax;
*
* $Header$
*/

class FullShowerProfileParamsManager{
 public:
  double energy;
  double alpha;
  double lnalpha;
  double slnalpha;
  double tmax;
  double lntmax;
  double slntmax;
  double rho;
  double relerr;
  double beta;
 private:
  double LNA0;
  double LNA1;
  double SLNA0;
  double SLNA1;
  double LNT0;
  double SLNT0;
  double SLNT1;
  double RHO0;
  double RHO1;
  double RHO2;
  double RELERR0;
  double RELERR1;
  double A0;
  double A1;
  double T0;
  double T1;

 public:
  FullShowerProfileParamsManager();
  virtual ~FullShowerProfileParamsManager();
  void Reset();
  bool Fill(double energy_input);
  double GetChi2Contribution(double *par);
};

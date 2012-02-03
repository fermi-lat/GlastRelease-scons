/**   
*
* THIS CODE IS STILL NOT FULLY DOCUMENTED !!!!!!!!!!!!!!!!!!!!
* DESCRIPTION COMMENTS WILL BE ADDED SOON.
* Philippe Bruel
*
*/

/**   
* @class NewFullShowerProfileParamsManager
* @author Philippe Bruel
*
* Class handling the energy parameterization of the shower profile parameters alpha and tmax;
*
* $Header$
*/

class NewFullShowerProfileParamsManager{
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
  double relerrparl;
  double relerrparr;
  double beta;
  double logenergy;
  double xcor;
  double ycor;
  int optmodprof;
  double modparam0;
  double modparam1;
  double modparam2;
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

  double THETACOR;
  double CTHETACOR;
  double STHETACOR;
  
  double XCORMEAN;
  double YCORMEAN;
  double XCORSIGMA;
  double YCORSIGMA;
  double XCORMPAR[3];
  double XCORSPAR[4];
  double YCORMPAR[3];
  double YCORSPAR[4];
  double XYERRPAR[4];
  double XYERRPARL[2];
  double XYERRPARR[3];
  double MODPROF00;
  double MODPROF01;
  double MODPROF02;
  double MODPROF20;
  double MODPROF21;
  double MODPROF22;
  double MODPROF23;
  double MODPROF24;

 public:
  NewFullShowerProfileParamsManager(int optmodprofin);
  virtual ~NewFullShowerProfileParamsManager();
  void Reset();
  bool Fill(double *par);
  double GetChi2Contribution(double *par);
};

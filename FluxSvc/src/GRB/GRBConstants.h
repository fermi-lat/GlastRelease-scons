#include <iostream.h>
#include <fstream.h>
#include "FluxSvc/mainpage.h"

#ifndef GRBCONSTANTS_HH
#define GRBCONSTANTS_HH 1

namespace cst
{
  /// Universal constants :
  const double mpc2 = 938.2;
  const double erg2MeV =624151.0;
  const double mpc2cm=3.0857e+24;
  //1Gauss = G2Mev * sqrt(Mev)
  const double G2MeV=815.78;
  const double st = 6.65225e-25;
  const double mec2 = 0.510999;
  const double c = 2.98e+10;
  const double c2 = c*c;
  const double hplanck=4.13567e-15; //eV*sec  
  const double pi = 3.1415926535897932385;
  const double Hubble=6.5e+1; 
  const double wzel=0.0;
  const double gamma0 = 100.;
  const double dgamma = 200.;
  const double csi = 0.1;
  const double alphae = .33;
  const double alphab = .33;
  const double p = 2.5;
  const double viscosity = 0.;
  /// Internal Parameters
  const double enmax = 1.0e+12;
  const double enph = 1.0e+7; /* Minimum photon energy detectable by GLAST */
  const double enmin = 1.0e+3;
  const double dt1=10000.;
  const int nstep=100;
  const int enstep=100;
  const double ch1L=2.0e+3;
  const double ch1H=5.0e+3;
  const double ch2L=5.0e+4;
  const double ch2H=3.0e+5;
  const double ch3L=2.5e+5;
  const double ch3H=1.8e+6;
  const double ch4L=2.0e+7;
  const double ch4H=3.0e+11;
  const float flagIC=0.1; // flag =[0,1], if ==0, No inverse compton;
}

class GRBConstants 
{ 
  
 public:
  GRBConstants(char filen);
  ~GRBConstants() { }
  
  void ReadParam(char filen);
  //  void GRBUniverse(double redshift);
  
  inline int Nshell() {return nshell;}
  inline void setNshell(int value=10){nshell=value;}
  
  inline double Redshift() {return redshift;}
  inline void setRedshift(double value){redshift=value;}
  
  inline double Etot(){return etot;}
  inline void setEtot(double value){etot=value;}
  
  inline double R0(){return r0;}
  inline void setR0(double value){r0=value;}
  
  inline double T0(){return t0;}
  inline void setT0(double value){t0=value;}

  int nshell;
  double redshift;
  double etot;
  double r0;
  double t0;

};

#endif



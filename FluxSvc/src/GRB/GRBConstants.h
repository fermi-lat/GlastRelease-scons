/**
 * GRBConstant: This file contains all the parameters 
 *              and all the constants needed for the GRB simulation program
 *              Some of them could be changed to run different GRB 
 *              descriptions. 
 * /author Nicola Omodei nicola.omodei@pi.infn.it 
 * /author Johann Coen Tanugi johann.cohen@pi.infn.it
 *
 * \b revision:
 *   02/15/2002
 *
 */
#include <iostream>
#include "FluxSvc/mainpage.h"
#ifndef GRBCONSTANTS_H
#define GRBCONSTANTS_H 1

class GRBConstants {
public:
  GRBConstants() { }
  ~GRBConstants() { }
  //  inline double Pi() { return acos(-1);}
  /// Universal constants :
  static const double mpc2 = 938.2;
  static const double erg2MeV =624151.0;
  static const double mpc2cm=3.0857e+24;
  //1Gauss = G2Mev * sqrt(Mev)
  static const double G2MeV=815.78;
  static const double st = 6.65225e-25;
  static const double mec2 = 0.510999;
  static const double c = 2.98e+10;
  static const double c2 = c*c;
  static const double hplanck=4.13567e-15; //eV*sec  
  static const double pi = 3.1415926535897932385;
  static const double Hubble=6.5e+1; 
  static const double wzel=0.0;
  /// Model Description:
  static const double red=0.1;
  static const int nshell=50;
  static const double etot=1.0e+55;
  static const double r0 = 1e+9;      //initial separation
  static const double dr0 = 0.1;  //r0*deltar0=initial, thikness
  static const double gamma0 = 100.;
  static const double dgamma = 200.;
  static const double csi = 0.1;
  static const double alphae = .33;
  static const double alphab = .33;
  static const double p = 2.5;
  static const double viscosity = 0.;
  /// Internal Parameters
  static const double enmax = 1.0e+12;
  static const double enph = 1.0e+7;
  static const double enmin = 1.0e+3;
  static const double dt1=10000.;

  static const int nstep=100;
  static const int enstep=100;

  static const double ch1L=2.0e+3;
  static const double ch1H=5.0e+3;
  static const double ch2L=5.0e+4;
  static const double ch2H=3.0e+5;
  static const double ch3L=2.5e+5;
  static const double ch3H=1.8e+6;
  static const double ch4L=2.0e+7;
  static const double ch4H=3.0e+11;
  static const float flagIC=0.1; // flag =[0,1], if ==0, No inverse compton;
  };
#endif








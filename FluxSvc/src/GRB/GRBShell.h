/**
 * GRBShell: Class that describes a Shell
 *
 * /author Nicola Omodei nicola.omodei@pi.infn.it 
 * /author Johann Coen Tanugi johann.cohen@pi.infn.it
 *
 * \b revision:
 *   02/15/2002
 *
 */
#include "FluxSvc/mainpage.h"

#include "GRBConstants.h"

#ifndef GRBSHELL_H
#define GRBSHELL_H 1

class GRBShell {
  
  static const GRBConstants cst;
  
 public:
  GRBShell(double);
  ~GRBShell() { }

  //Accessors:
  inline double Mass() {return _mass;}
  inline double Gamma() {return _gamma;}
  inline double Thickness() {return _thickness;}
  inline double Radius() {return _radius;}
  inline double VolCom() {return 4.0*cst.pi*_thickness*_radius*_radius*_gamma;}

  //Set functions
  inline void setMass(double value) {_mass = value;}
  inline void setGamma(double value) {_gamma = value;}
  inline void setThickness(double value) {_thickness = value;}
  inline void setRadius(double value) {_radius = value;}
  
  //Higher  level functions:
  double    generateGamma(double gamma0,double dgamma);
  double    beta(const double gamma);
  void      evolve(double time);


 private:
  //Data Members:
  double _mass;
  double _gamma;
  double _thickness;
  double _radius;
  int    _status;
};

#endif


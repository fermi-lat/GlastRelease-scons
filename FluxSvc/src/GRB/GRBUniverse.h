/**
 * GRBUniverse:Access to some cosmological quantities
 *
 * /author Nicola Omodei nicola.omodei@pi.infn.it 
 * /author Johann Coen Tanugi johann.cohen@pi.infn.it
 *
 * \b revision:
 *   02/15/2002
 *
 */
#include "FluxSvc/mainpage.h"
#ifndef GRBUNIVERSE_HH
#define GRBUNIVERSE_HH 1

class GRBUniverse {
  GRBConstants cst;
public:
  GRBUniverse();
  ~GRBUniverse() { }
  inline double qo() {return _qo;}
  inline double Dist(){return _Dist;}
  inline double Area(){return _Area;}
  double _qo;
  double _Dist;
  double _Area;
  double _Hubble;
};
#endif

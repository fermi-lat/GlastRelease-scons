#ifndef GDDMPOSY_H
#define GDDMPOSY_H
#include "detModel/Sections/GDDanyPosition.h"

/**
 * @author D.Favretto & R.Giannitrapani 
 */
class GDDmposY : GDDanyPosition {
public:

  GDDmposY( int pncopy,double pdy,positiontype ptypePosition);
  /**
   *
   */ 
  void setY0(double py0);
  /**
   *
   */ 
  void setZX(double pz,double px);
 private:    
  int ncopy;
  double y0; ///default 0
  double dy;
  double z; ///default 0
  double x; ///default 0
};
#endif //GDDMPOSY_H


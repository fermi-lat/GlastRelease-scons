#ifndef GDDMPOSX_H
#define GDDMPOSX_H
#include "detModel/Sections/GDDanyPosition.h"

/**
 * @author D.Favretto & R.Giannitrapani 
 */
class GDDmposX : GDDanyPosition {
public:

  GDDmposX(int pncopy,double pdx,positiontype ptypePosition);
  /**
   *
   */ 
  void setX0(double px0);
  /**
   *
   */ 
  void setYZ(double pyz[]);
  
 private:    
  int ncopy;
  double x0; ///default 0
  double dx;
  double yz[2]; ///default "0 0"
};
#endif //GDDMPOSX_H

#ifndef GDDMPOSPHI_H
#define GDDMPOSPHI_H
#include "detModel/Sections/GDDanyPosition.h"

/**
 * @author D.Favretto & R.Giannitrapani 
 */
class GDDmposPhi : GDDanyPosition {
public:

  GDDmposPhi(int pncopy,double pdPhi,positiontype ptypePosition);
  /**
   *
   */ 
  void setRot(bool pimpliedRot);
  /**
   *
   */ 
  void setPhi0(double pPhi0);
  /**
   *
   */ 
  void setRZ(double pRZ[]);
  
 private:    
  int ncopy;
  double Phi0; ///defaut 0
  double dPhi;
  double RZ[2]; ///defaul "0 0"
  bool impliedRot; ///default true
};
#endif //GDDMPOSPHI_H


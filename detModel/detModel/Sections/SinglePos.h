#ifndef SINGLEPOS_H
#define SINGLEPOS_H
#include "detModel/Sections/Position.h"

namespace detModel{

  /**
   * @author R.Giannitrapani
   */
  class SinglePos :public Position {
  public:
    SinglePos():Position(),xRot(0),yRot(0),zRot(0),S(0){;}

    /**
     *  This methods sets and gets the xRot,yRot,zRot,S
     */ 
    double getXRot()const{ return xRot; }
    void setXRot(double XRot){xRot = XRot; }
    double getYRot()const{ return yRot; }
    void setYRot(double YRot){yRot = YRot; }
    double getZRot()const{ return zRot; }
    void setZRot(double ZRot){zRot = ZRot; }
    double getS()const{ return S; }
    void setS(double s){ S = s; }

  private:  
    double xRot;
    double yRot;
    double zRot;
    double S;
  };
}
#endif //SINGLEPOS_H

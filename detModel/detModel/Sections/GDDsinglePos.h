#ifndef GDDSINGLEPOS_H
#define GDDSINGLEPOS_H
#include "detModel/Sections/GDDposition.h"

/**
 * @author R.Giannitrapani
 */
class GDDsinglePos :public GDDposition {
 public:
  GDDsinglePos():GDDposition(),xRot(0),yRot(0),zRot(0),S(0){;}

  /**
   *  This methods sets and gets the xRot,yRot,zRot,S
   */ 
  double getXRot()const{ return XRot; }
  void setXRot(double XRot){xRot = XRot; }
  double getYRot()const{ return YRot; }
  void setYRot(double YRot){yRot = YRot; }
  double getZRot()const{ return yRot; }
  void setZRot(double ZRot){zRot = ZRot; }
  double getS()const{ return S; }
  void setS(double s){ S = s; }

private:  
  double xRot;
  double yRot;
  double zRot;
  double S;
};
#endif //GDDSINGLEPOS_H

#ifndef GDDSTACKEDPOS_H
#define GDDSTACKEDPOS_H
#include "detModel/Sections/GDDposition.h"

/**
 * @author R.Giannitrapani
 */
class GDDstackedPos :public GDDposition {
 public:
  enum axisDir{xDir, yDir, zDir};

 public:
  GDDstackedPos(axisDir dir):
    GDDposition(),aDir(dir),
    dx(0),dy(0),dz(0),
    rotation(0),shift(0),gap(0){;}

  GDDstackedPos():
    GDDposition(),aDir(xDir),
    dx(0),dy(0),dz(0),
    rotation(0),shift(0),gap(0){;}
  
  axisDir getAxisDir(){return aDir;};
  void setAxisDir(axisDir dir){aDir = dir;};

  void setDx(double pdx){dx=pdx;}
  void setDy(double pdy){dy=pdy;}
  void setDz(double pdz){dz=pdz;}
  double getDx(){return dx;}
  double getDy(){return dy;}
  double getDz(){return dz;}

  void setRotation(double protation){rotation=protation;}
  double getRotation(){return rotation;}

  void setShift(double pshift){shift=pshift;}
  void setGap(double pgap){gap=pgap;}
  double getShift()const{return shift;}
  double getGap()const{return gap;}


 private:
  /// The direction of the stacked positioning
  axisDir aDir;

  double dx; 
  double dy; 
  double dz; 
  double rotation; 
  double shift; 
  double gap;  
};
#endif //GDDSTACKEDPOS_H

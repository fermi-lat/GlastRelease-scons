#ifndef GDDPOSXYZ_H
#define GDDPOSXYZ_H
#include "detModel/Sections/GDDsinglePos.h"
#include "detModel/Management/GDDsectionsVisitor.h"

/**
 * @author D.Favretto & R.Giannitrapani
 */

class GDDposXYZ :public GDDsinglePos {
 public:
  /// A standard constructor
  GDDposXYZ():GDDsinglePos(),x(0),y(0),z(0){;}
  /// A constructor with the coordinate of the positioned volume
  GDDposXYZ(double x, double y, double z):GDDsinglePos(),
    x(x),y(y),z(z){;}

  void Accept(GDDsectionsVisitor* v);
  virtual void AcceptNotRec(GDDsectionsVisitor* v){v->visitPosXYZ(this);};

  void setX(double pX){x=pX;}
  void setY(double pY){y=pY;}
  void setZ(double pZ){z=pZ;}
  double getX()const{return x;}
  double getY()const{return y;}
  double getZ()const{return z;}

  virtual void buildBB();
 private:    
  double x; 
  double y; 
  double z; 
  
};
#endif //GDDPOSXYZ_H




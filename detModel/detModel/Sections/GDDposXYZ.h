#ifndef GDDPOSXYZ_H
#define GDDPOSXYZ_H
#include "detModel/Sections/GDDanySinglePosition.h"
#include "detModel/Management/GDDsectionsVisitor.h"

/**
 * @author D.Favretto & R.Giannitrapani
 */

class GDDposXYZ :public GDDanySinglePosition {
public:
  GDDposXYZ(positiontype ptypePosition):GDDanySinglePosition(ptypePosition){;}
  /**
   *
   */ 
  void Accept(GDDsectionsVisitor* v);
  /**
   *
   */ 
  virtual void AcceptNotRec(GDDsectionsVisitor* v){v->visitPosXYZ(this);};
  /**
   *
   */ 
  void setX(double pX){x=pX;}
  /**
   *
   */ 
  void setY(double pY){y=pY;}
  /**
   *
   */ 
  void setZ(double pZ){z=pZ;}
  /**
   *
   */ 
  virtual double getBBX();
  /**
   *
   */ 
  virtual double getBBY();
  /**
   *
   */ 
  virtual double getBBZ();
  /**
   *
   */ 
  double getX()const{return x;}
  /**
   *
   */ 
  double getY()const{return y;}
  /**
   *
   */ 
  double getZ()const{return z;}
  
 private:    
  double x; ///default 0  NMTOKEN
  double y; ///default 0  NMTOKEN
  double z; ///default 0  NMTOKEN
  
};
#endif //GDDPOSXYZ_H




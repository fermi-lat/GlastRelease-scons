#ifndef GDDAXISPOS_H
#define GDDAXISPOS_H
#include "detModel/Sections/GDDanyRelativePosition.h"
#include "detModel/Management/GDDsectionsVisitor.h"

/**
 * @author D.Favretto & R.Giannitrapani
 */
class GDDaxisPos :public GDDanyRelativePosition {
 public:

  GDDaxisPos():
    GDDanyRelativePosition(AxisPos),shift(0),gap(0){}

  virtual ~GDDaxisPos(){};
  /**
   *
   */   
  virtual void Accept(GDDsectionsVisitor* v);
  /**
   *
   */ 
  virtual void AcceptNotRec(GDDsectionsVisitor* v){v->visitAxisPos(this);};
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
  void setShift(double pshift){shift=pshift;}
  /**
   *
   */ 
  void setGap(double pgap){gap=pgap;}
  /**
   *
   */ 
  double getShift()const{return shift;}
  /**
   *
   */ 
  double getGap()const{return gap;}
  
 private:
  double shift; ///default 0
  double gap;  ///default 0
};
#endif //GDDAXISPOS_H

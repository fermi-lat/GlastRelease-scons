#ifndef GDDAXISMPOS_H
#define GDDAXISMPOS_H
#include "detModel/Sections/GDDanyRelativePosition.h"
#include "detModel/Management/GDDsectionsVisitor.h"

/**
 * @author D.Favretto & R.Giannitrapani 
 */
class GDDaxisMPos : public GDDanyRelativePosition {
public:

  GDDaxisMPos():
    GDDanyRelativePosition(AxisMpos),shift0(0),gap0(0),ncopy(0){};
  /**
   *
   */ 
  virtual void Accept(GDDsectionsVisitor* v);  
  /**
   *
   */ 
  virtual void AcceptNotRec(GDDsectionsVisitor* v){v->visitAxisMPos(this);};
  /**
   *
   */ 
  void setShift0(double pshift0){shift0=pshift0;};
  /**
   *
   */ 
  void setGap0(double pgap0){gap0 = pgap0;};
  /**
   *
   */ 
  void setNcopy(int pncopy){ncopy = pncopy;};
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
  double getShift0()const{return shift0;}
  /**
   *
   */ 
  double getGap0()const{return gap0;}
  /**
   *
   */ 
  int getNcopy()const{return ncopy;}

private:    
    double shift0; ///default 0
    double gap0; ///default 0
    int ncopy; ///default 0

};
#endif //GDDAXISMPOS_H








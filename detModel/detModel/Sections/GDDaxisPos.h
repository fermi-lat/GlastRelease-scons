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
    GDDanyRelativePosition(AxisPos){};
  
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

 private:
  
  // This is the displacement along the axis of the stack of the
  // center of mass of the positioned volume; the visitors can use
  // this quantity to pose in the right way the volume.
  double disp;

};
#endif //GDDAXISPOS_H


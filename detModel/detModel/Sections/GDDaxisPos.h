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
    GDDanyRelativePosition(AxisPos){}

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


};
#endif //GDDAXISPOS_H

#ifndef GDDAXISPOS_H
#define GDDAXISPOS_H
#include "detModel/Sections/GDDstackedPos.h"
#include "detModel/Management/GDDsectionsVisitor.h"

/**
 * @author D.Favretto & R.Giannitrapani
 */
class GDDaxisPos :public GDDstackedPos {
 public:

  GDDaxisPos(axisDir dir):GDDstackedPos(dir),disp(0){};
  GDDaxisPos():disp(0){};
  
  virtual ~GDDaxisPos(){};
  virtual void Accept(GDDsectionsVisitor* v);
  virtual void AcceptNotRec(GDDsectionsVisitor* v){v->visitAxisPos(this);};

  void setDisp(double d){disp = d;};
  double getDisp(){return disp;};

  virtual void buildBB();

 private:
  
  // This is the displacement along the axis of the stack of the
  // center of mass of the positioned volume; the visitors can use
  // this quantity to pose in the right way the volume.
  double disp;

};
#endif //GDDAXISPOS_H


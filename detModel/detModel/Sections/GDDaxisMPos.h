#ifndef GDDAXISMPOS_H
#define GDDAXISMPOS_H
#include <vector>
#include "detModel/Sections/GDDstackedPos.h"
#include "detModel/Management/GDDsectionsVisitor.h"

/**
 * @author D.Favretto & R.Giannitrapani 
 */
class GDDaxisMPos : public GDDstackedPos {
public:

  GDDaxisMPos(axisDir dir):GDDstackedPos(dir),shift0(0),gap0(0),ncopy(0){};
  GDDaxisMPos():GDDstackedPos(),shift0(0),gap0(0),ncopy(0){};

  virtual void Accept(GDDsectionsVisitor* v);  
  virtual void AcceptNotRec(GDDsectionsVisitor* v){v->visitAxisMPos(this);};

  void setShift0(double pshift0){shift0=pshift0;};
  void setGap0(double pgap0){gap0 = pgap0;};
  void setNcopy(int pncopy){ncopy = pncopy;};
  double getShift0(){return shift0;}
  double getGap0(){return gap0;}
  int getNcopy(){return ncopy;}

  virtual void buildBB();
  double getDispCM(){return dispCM;};
  void setDispCM(double d){dispCM = d;};
  void addDisp(double d){disp.push_back(d);};
  void setDisp(int i, double d){disp[i] = d;}
  double getDisp(int i){return disp[i];};
  void translateDisp(double d);
private:    
  double shift0; ///default 0
  double gap0; ///default 0
  int ncopy; ///default 0

  double dispCM;

  std::vector <double> disp;
};
#endif //GDDAXISMPOS_H








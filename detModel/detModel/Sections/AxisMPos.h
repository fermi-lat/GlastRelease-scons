#ifndef AXISMPOS_H
#define AXISMPOS_H
#include <vector>
#include "detModel/Sections/StackedPos.h"
#include "detModel/Management/SectionsVisitor.h"

namespace detModel{

  /**
   * This is a multiple axis position
   * @author D.Favretto & R.Giannitrapani 
   */
  class AxisMPos : public StackedPos {
  public:

  
    AxisMPos(Stack::axisDir dir):StackedPos(dir),shift0(0),gap0(0),ncopy(0){};
    AxisMPos():StackedPos(),shift0(0),gap0(0),ncopy(0){};

    virtual void Accept(SectionsVisitor* v);  
    virtual void AcceptNotRec(SectionsVisitor* v){v->visitAxisMPos(this);};

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
}
#endif //AXISMPOS_H








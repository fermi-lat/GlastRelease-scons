#ifndef STACK_H
#define STACK_H
#include <vector>
#include "detModel/Sections/Ensemble.h"

namespace detModel{

  /**
   * @author D.Favretto & R.Giannitrapani
   */

  class Stack :public Ensemble {
  public:
    enum stackOrigin{atStart,atCentre};
    enum axisDir{xDir, yDir, zDir};
  
  public:
    Stack(std::string pName,axisDir dir):Ensemble(pName),aDir(dir),origin(atStart){;}
    Stack(axisDir dir):Ensemble(),aDir(dir),origin(atStart){;}
    virtual ~Stack(){;}

    void setOrigin(stackOrigin porigin){origin=porigin;}
    stackOrigin getOrigin()const{return origin;}
  
    void setAxisDir(axisDir dir){aDir = dir;};
    axisDir getAxisDir()const{return aDir;}

    virtual void buildBB();

    /// This method build the displacement list of its positioning
    void buildDisp();

    /// Stack uses this Template Method to put a check on position to be added
    bool checkPosition(Position* p);

  private:
    axisDir aDir;
    stackOrigin origin; ///default atStart
  };
}
#endif //STACK_H

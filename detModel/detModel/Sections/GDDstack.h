#ifndef GDDSTACK_H
#define GDDSTACK_H
#include <vector>

#include "GDDensemble.h"
#include "GDDensemble.h"

/**
 * @author D.Favretto & R.Giannitrapani
 */

class GDDstack :public GDDensemble {
 public:
  enum stackOrigin{atStart,atCentre};
  enum axisDir{xDir, yDir, zDir};
  
 public:
  GDDstack(std::string pName,axisDir dir):GDDensemble(pName),aDir(dir),origin(atStart){;}
  GDDstack(axisDir dir):GDDensemble(),aDir(dir),origin(atStart){;}
  virtual ~GDDstack(){;}

  void setOrigin(stackOrigin porigin){origin=porigin;}
  stackOrigin getOrigin()const{return origin;}
  
  void setAxisDir(axisDir dir){aDir = dir;};
  axisDir getAxisDir()const{return aDir;}

  virtual void buildBB();

  /// This method build the displacement list of its positioning
  void buildDisp();

  /// Stack uses this Template Method to put a check on position to be added
  bool checkPosition(GDDposition* p);

 private:
  axisDir aDir;
  stackOrigin origin; ///default atStart
};
#endif //GDDSTACK_H









#ifndef GDDSTACK_H
#define GDDSTACK_H
#include <vector>

#include "detModel/Sections/GDDensamble.h"

/**
 * @author D.Favretto & R.Giannitrapani
 */

class GDDstack :public GDDensamble {
 public:
  enum stackOrigin{atStart,atCentre};
  enum axisDir{xDir, yDir, zDir};
  
 public:
  GDDstack(std::string pName,axisDir dir):GDDensamble(pName),aDir(dir),origin(atStart){;}
  GDDstack(axisDir dir):GDDensamble(),aDir(dir),origin(atStart){;}
  virtual ~GDDstack(){;}

  void setOrigin(stackOrigin porigin){origin=porigin;}
  stackOrigin getOrigin()const{return origin;}
  
  void setAxisDir(axisDir dir){aDir = dir;};
  axisDir getAxisDir()const{return aDir;}

  virtual void buildBB();

  /// This method build the displacement list of its positioning
  void buildDisp();

 private:
  axisDir aDir;
  stackOrigin origin; ///default atStart
};
#endif //GDDSTACK_H









#ifndef GDDBOUNDABLE_H
#define GDDBOUNDABLE_H

#include "detModel/Sections/GDDboundingBox.h"

/**
 * This class represent a generic object that can have
 * a bounding box (i.e. a volume or a positioning); this
 * simplify the abstractions in detModel
 * @author R.Giannitrapani
 */
class GDDboundable {
 public:
  GDDboundable(){bBox = new GDDboundingBox();}
  ~GDDboundable(){delete bBox;};
  /**
   * This abstract method build the bounding box.
   * It has to be concretely implemented in sublcass
   */
  virtual void buildBB() = 0;

  /// This method return the pointer to the bounding box
  GDDboundingBox* getBBox(){return bBox;};

 private:
  /// The bounding box 
  GDDboundingBox* bBox;
};
#endif //GDDBOUNDABLE_H

#ifndef BOUNDABLE_H
#define BOUNDABLE_H
#include "detModel/Sections/BoundingBox.h"

namespace detModel{

  /**
   * This class represent a generic object that can have
   * a bounding box (i.e. a volume or a positioning); this
   * simplify the abstractions in detModel
   * @author R.Giannitrapani
   */
  class Boundable {
  public:
    Boundable(){bBox = new BoundingBox();}
    virtual ~Boundable(){delete bBox;};
    /**
     * This abstract method build the bounding box.
     * It has to be concretely implemented in sublcass
     */
    virtual void buildBB() = 0;
    
    /// This method return the pointer to the bounding box
    BoundingBox* getBBox(){return bBox;};
    
  private:
    /// The bounding box 
    BoundingBox* bBox;
  };
}
#endif //BOUNDABLE_H

#ifndef UNITS_H
#define UNITS_H
#include <string>

namespace detModel{

  /**
   * This class represents lengths and angles units
   * for shapes and positioned volumes
   * @author R.Giannitrapani
   */
  class Units {
  public:

    /// Enumerative for the lengths units
    enum unitLengthType{mm,cm,m};
    /// Enumerative for the angles units
    enum unitAngleType{deg,mrad};
  

  public:
    /**
     * This method set the units for the lenghts of the shape
     */  
    void setUnitLength(unitLengthType u){unitLength = u;};
    void setUnitLength(std::string s);
    /**
     * This method set the units for the angles of the shape
     */  
    void setUnitAngle(unitAngleType u){unitAngle = u;};
    void setUnitAngle(std::string s);


    /**
     * This method get the units for the lenghts of the shape
     */  
    unitLengthType getUnitLength(){return unitLength;};
    /**
     * This method get the units for the angles of the shape
     */  
    unitAngleType getUnitAngle(){return unitAngle;};

  private:
    unitLengthType unitLength; 
    unitAngleType unitAngle; 
  };
}
#endif //UNITS_H

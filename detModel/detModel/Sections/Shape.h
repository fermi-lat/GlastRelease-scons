#ifndef SHAPE_H
#define SHAPE_H
#include "detModel/Sections/Volume.h"
#include "detModel/Sections/Units.h"

namespace detModel{

  /**
   * This class contains all the common properties and 
   * functionalities of shapes. Concrete shapes must 
   * sublcass this one.
   * @author D.Favretto & R.Giannitrapani 
   */
  class Shape : public Volume {
  public:
    /** Some constructors; since this is an abstract class, these
     * constructors are useful only in its concrete subclasses
     */
    Shape(std::string pname):Volume(pname),sensitive(0){units = new Units;}
    Shape():Volume(),sensitive(0){units = new Units;}

    virtual ~Shape(){delete units;};

    /**
     * With this method it is possible to set the shape
     * sensitive for simulations purpouses
     */
    void setSensitive(bool psensitive){sensitive = psensitive;};
    bool getSensitive(){return sensitive;};
    /**
     * This method set the material (a string) of the shape
     */  
    void setMaterial(std::string pmaterial){material = pmaterial;};
    /**
     * This method gives back the material (a string) of the shape
     */  
    std::string getMaterial(){return material;};

    Units* getUnits(){return units;};

  private:
    std::string material;
    bool sensitive;          
    Units* units;
  };
}
#endif //SHAPE_H

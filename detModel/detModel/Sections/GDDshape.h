#ifndef GDDSHAPE_H
#define GDDSHAPE_H
#include "detModel/Sections/GDDvolume.h"
#include "detModel/Sections/GDDunits.h"

/**
 * This class contains all the common properties and 
 * functionalities of shapes. Concrete shapes must 
 * sublcass this one.
 * @author D.Favretto & R.Giannitrapani 
 */
class GDDshape : public GDDvolume {
 public:
  /** Some constructors; since this is an abstract class, these
   * constructors are useful only in its concrete subclasses
   */
  GDDshape(std::string pname):GDDvolume(pname),sensitive(0){units = new GDDunits;}
  GDDshape():GDDvolume(),sensitive(0){units = new GDDunits;}

  virtual ~GDDshape(){delete units;};

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

  GDDunits* getUnits(){return units;};

 private:
  std::string material;
  bool sensitive;          
  GDDunits* units;
};
#endif //GDDSHAPE_H









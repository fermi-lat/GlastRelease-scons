#ifndef GDDSHAPE_H
#define GDDSHAPE_H
#include "detModel/Sections/GDDvolume.h"

class GDDsectionsVisitor;

/// Enumerative for the identification of the type of shape 
enum shapeType{box, tube, cone};

/**
 * This class contains all the common properties and 
 * functionalities of shape. Concrete shapes must 
 * hinerith from it.
 * @author D.Favretto & R.Giannitrapani 
 */
class GDDshape : public GDDvolume {
 public:

  /** 
   * This is the standard Accept method for the visitor mechanism
   */
  virtual void Accept(GDDsectionsVisitor*){};

  /**
   * With this method it is possible to set the shape
   * sensitive for simulations purpouses
   */
  void setSensitive(bool psensitive){sensitive = psensitive;};

  bool getSensitive(){return sensitive;};

  /**
   * This method set the material (a string) of the shape
   */  
  void setMaterial(string pmaterial){material = pmaterial;};

  /**
   * This method gives back the material (a string) of the shape
   */  
  string getMaterial(){return material;};


  /**
   * This method set the units for the lenghts of the shape
   */  
  void setUnitLength(char* punitLength){

    if (punitLength== "mm")
      unitLength = GDDmm;
    else if (punitLength== "cm")
      unitLength = GDDcm;
    else unitLength = GDDm;
  };

  /**
   * This method set the units for the angles of the shape
   */  
  void setUnitAngle(char* punitAngle){
    if (punitAngle== "deg")
      unitAngle = GDDdeg;
    else unitAngle = GDDmrad;
  };

  /**
   * This method get the units for the lenghts of the shape
   */  
  unitLengthtype getUnitLength(){return unitLength;};

  /**
   * This method get the units for the angles of the shape
   */  
  unitAngletype getUnitAngle(){return unitAngle;};
  /**
   * This method gives the type of the shape
   */  
  shapeType getShapeType(){return stype;};

 protected:
  GDDshape(string pname, shapeType pstype):GDDvolume(pname,shape),stype(pstype){;}
  GDDshape(){};
  virtual ~GDDshape(){};

 private:
  shapeType stype;
  string material;
  bool sensitive;          ///default false
  unitLengthtype unitLength; ///dufault mm
  unitAngletype unitAngle; ///default deg
};
#endif //GDDSHAPE_H









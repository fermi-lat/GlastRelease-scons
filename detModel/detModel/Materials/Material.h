#ifndef MATERIAL_H
#define MATERIAL_H
#include <vector>
#include <string>
#include "detModel/Utilities/Color.h"

namespace detModel{

  class MaterialsVisitor;

  /**
   * This class represents a material for shapes of GLAST geometry;
   * it is an abstrat class.
   * @author M.Frailis & R.Giannitrapani 
   */
  class Material{
  public:
    
    /// A standard constractor
    Material():name(""),density(0.0){}  

    /// A standard constractor
    Material(std::string n, double d):name(n),density(d) {;}  
    
    /// A standard destructor
    virtual ~Material(){}

    /// Methods to set and get the name
    void setName(std::string n){name=n;}  
    std::string  getName(){return name;}

    /// Methods to set and get the density
    void setDensity(double d){density=d;}  
    double getDensity(){return density;}

    /// This is the standard Accept method for the visitor mechanism
    virtual void Accept(MaterialsVisitor*) = 0;
    /// This method is the non recursive Accept for the visitor pattern
    virtual void AcceptNotRec(MaterialsVisitor*) = 0;
    
  private:
    /// The material's name
    std::string name;
    /// The material's density
    double density;
  };
}
#endif //MATERIAL_H








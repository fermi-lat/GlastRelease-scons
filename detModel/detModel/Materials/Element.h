#ifndef ELEMENT_H
#define ELEMENT_H
#include <vector>
#include <string>
#include "detModel/Materials/Material.h"

namespace detModel{

  class MateriasVisitor;

  /**
   * This class represents an element used in materials;
   * @author M.Frailis & R.Giannitrapani 
   */
  class Element: public Material{
  public:

    /// A standard constractor
    Element(): Material() {symbol=""; z=0; aweight=0.0;}  
    
    /// A standard destructor
    virtual ~Element(){;}
    
    /// Methods to set and get the symbol
    void setSymbol(std::string s){symbol=s;}  
    std::string getSymbol(){return symbol;}

    /// Methods to set and get the z value
    void setZ(int t){z=t;}  
    int getZ(){return z;}

    /// Methods to set and get the average weight
    void setAweight(double a){aweight=a;}  
    double getAweight(){return aweight;}
    
    /// This is the standard Accept method for the visitor mechanism
    virtual void Accept(MaterialsVisitor*);
    /// This method is the non recursive Accept for the visitor pattern
    virtual void AcceptNotRec(MaterialsVisitor*);
    
  private:
     
    /// The element's symbol
    std::string symbol;
    /// The element's z value
    int z;
    /// The element's average weight
    double aweight;

  };
}
#endif //ELEMENT_H





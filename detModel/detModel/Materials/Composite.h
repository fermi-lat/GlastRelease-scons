#ifndef COMPOSITE_H
#define COMPOSITE_H
#include <vector>
#include <string>
#include "detModel/Utilities/Purge.h"
#include "detModel/Utilities/Global.h"
#include "detModel/Materials/Material.h"

namespace detModel{
  
  class MaterialsVisitor;
  
  /**
   * This class represents a composite material for shapes of GLAST geometry;
   * it is built by adding element or other composite with the aid of the
   * abstract class Material. A composite is composed by a list of material
   * specifing either the number of atoms of each material or its fraction.
   * @author M.Frailis & R.Giannitrapani 
   */
  class Composite: public Material{
  public:
    /// A constructor with the number of components
    Composite(unsigned int n): Material(),nComponents(n) {}  
    
    /// A standard destructor
    virtual ~Composite(){;}

    /// Return the number of components that this material should have
    unsigned int getNComponents(){return nComponents;}

    /// These methods return the vectors of components, they atoms number and their fraction 
    std::vector<Material*> getComponents()const{return components;}
    std::vector<int> getAtoms()const{return atoms;}
    std::vector<double> getFractions()const{return fractions;}


    /// This is the standard Accept method for the visitor mechanism
    virtual void Accept(MaterialsVisitor*);
    /// This method is the non recursive Accept for the visitor pattern
    virtual void AcceptNotRec(MaterialsVisitor*);


    /** 
     *  This methods return the atom number and the fraction of a given material
     *  If the name does not correspond to nothing, 0 is returned
     */
    unsigned int getAtom(std::string mName);
    double getFraction(std::string mName);
    
    /// These methods return true if the components are specified by fractions or atoms
    bool isFractions(){return fractions.size();}
    bool isAtoms(){return atoms.size();}

    /// This overloaded method add a material to the list of components
    void addComponent(Material*, unsigned int);
    void addComponent(Material*, double);

  private:
    /// The number of components
    unsigned int nComponents;
    /// The componentes vector
    std::vector <Material*> components;
    /// The vector with number of atoms for components
    std::vector <int> atoms;
    /// The vector with fractions of the components
    std::vector <double> fractions;    
  };
}
#endif //COMPOSITE_H




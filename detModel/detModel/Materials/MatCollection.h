#ifndef MATCOLLECTION_H
#define MATCOLLECTION_H
#include <vector>
#include <map>
#include <string>
#include "detModel/Utilities/Purge.h"
#include "detModel/Materials/Material.h"

namespace detModel{

  class MaterialsVisitor;
  class Material;
  class Color;

  /**
   * This class represents a collection of materials;
   * it is suited for use with the gdd.dtd materials element.
   * @author M.Frailis & R.Giannitrapani 
   */
  class MatCollection{
  public:

    /// A standard constructor
    MatCollection():version(""),date(""),author(""){;}  
    /// The destructor call the purge for the materials vector
    virtual ~MatCollection(){}

    /// Methods to set and get the version
    void setVersion(std::string v){version=v;}  
    std::string getVersion(){return version;}

    /// Methods to set and get the date
    void setDate(std::string d){date=d;}  
    std::string getDate(){return date;}

    /// Methods to set and get the author
    void setAuthor(std::string a){author=a ;}  
    std::string getAuthor(){return author;}

    
    /// Method to add a material name to the collection
    void addMaterialName(std::string);
    /// Method to add a material to the collection
    void addMaterial(Material*);
    /// Method to get the total number of materials in the collection
    int getMaterialsNumber(){return materials.size();};
    /// Method to retrive a materials map
    std::map<std::string, Material* > getMaterials(){return materials;};
    /// Method to retrive a material colors map
    std::map<std::string, Color* > getMaterialColors(){return materialColors;};

    
    /// This is the standard Accept method for the visitor mechanism
    void Accept(MaterialsVisitor*);
    /// This method is the non recursive Accept for the visitor pattern
    void AcceptNotRec(MaterialsVisitor*);

    /// This method redistributes the colors of material not fixed
    void redistributeColor(){};
    /// This method assign automatically all the colors
    void generateColor();
    
    /// This method sets if the material color is fixed or not
    void setFixedMaterial(std::string, bool){}; 

    /// This methods set a material color directly
    void setMaterialColor(std::string, float, float, float, float);
    void setMaterialColor(std::string, Color*);

    /// This method gives back a material color
    Color* getMaterialColor(std::string);

    /// This methods set a material transparency directly
    void setMaterialTransparency(std::string,float);
    /// This methods set all the materials transaprency directly
    void setAllMaterialTransparency(float);

  private:    
    /// The version of the material collection
    std::string version;
    /// The date of the material collection
    std::string date;
    /// The author of the material collection
    std::string author;

    /// The map containing the colors for materials
    std::map <std::string, Color*> materialColors;
    /// The map containing the materials
    std::map <std::string, Material * > materials;
    /// The map containing the bool to say if a material color is fixed
    std::map <std::string, bool> fixedMaterial;    

  };
}
#endif //MATCOLLECTION_H








#ifndef GDD_H
#define GDD_H
#include <map>
#include <vector>
#include <string>
#include <list>
#include "detModel/Sections/Section.h"

namespace xmlUtil {
  class IdDict;
}


namespace detModel{
  
  using xmlUtil::IdDict;

  class Volume;
  class Constants;
  class Const;
  class ConstCategory;
  class Visitor;
  class Choice;
  class Shape;
  class MatCollection;
  class Material;
  class Color;


  /**
   * This is the main container of all the geometry.
   * Responsability of its creation and destruction is of the manager.
   * This class provides a series of maps and methods to query that 
   * maps. 
   * 
   * @author D.Favretto & R.Giannitrapani
   */
  class Gdd {    
  public:
    /** This is the destructor; it should be called only by the manager
     * destructor. It starts the destruction of all the objects created
     * by the builder
     */
    virtual ~Gdd();

    /// This method gives back the DTDversion
    std::string getDTDversion(){return DTDversion;};
    /// This method gives back the CVSid  
    std::string getCVSid(){return CVSid;};
    /// This method sets the DTDversion  
    void setDTDversion(std::string pdtd){DTDversion = pdtd;};
    /// This method sets the CVSid  
    void setCVSid(std::string pcvs){CVSid = pcvs;};

    ///This method sets the constants 
    void setConstants(Constants* pconstants){constants=pconstants;}
    ///This method gives back the constants
    Constants* getConstants(){return constants;}
    /// This method gives back a pointer to the sections vector
    std::vector <Section*>  getSections(){return sections;};  
    /// This method adds a section to the sections vector
    void addSection(Section* s){sections.push_back(s);};
    ///This method sets the materials collection
    void setMaterials(MatCollection* m){materials=m;}
    ///This method gives back the materials collection
    MatCollection* getMaterials(){return materials;}

    /// Keep pointer to id dictionary
    void setIdDictionary(IdDict* d) {dict=d;}
    /// Retrieve pointer to id dictionary
    IdDict* getIdDictionary() {return dict;}

    ///This method gives back the Gddvolume map  
    std::map < std::string, Volume * > getVolumesMap(){return volumeMap;};

    /// This is the recursive accept for the visitor pattern
    void Accept(Visitor* v);
    /// This is the non recursive accept for the visitor pattern
    void AcceptNotRec(Visitor* v);

    /// This method builds a global volumes map for all the sections 
    void buildVolumeMap();
    /// This method builds the global constants maps for all the constants
    void buildConstantsMap();
    
    /// This method builds the bounding boxes of all the boundable objects
    void buildBoundingBoxes(); 
    
    /// This method gives back the modes names vector
    std::vector <std::string> getModeNames(){return modeNames;};

    /**
     * This method search the volumes map with the name string and return 
     * a pointer to the Gddvolume if it exists, otherwise it returns a null 
     * pointer.
     */
    Volume * getVolumeByName(std::string vname);

    /**
     * This method search the volumes map with the name string and return 
     * a pointer to the Gddshape if it exists, otherwise it returns a null 
     * pointer. 
     */
    Shape * getShapeByName(std::string vname);
  
    /**
     * This method search the constants categories with the name string and return 
     * a pointer to the GddconstCategory if it exists, otherwise it returns a null 
     * pointer. 
     */
    ConstCategory* getConstCategoryByName(std::string vname);

    /**
     * This method search the const map with the name string and return
     * a pointer to the Gddconst object if it exists, otherwise it returns 0
     */
    Const* getConstByName(std::string cname);

    /**
     * This method search the materials map with the name string and return
     * a pointer to the Material object if it exists, otherwise it returns 0
     */
    Material* getMaterialByName(std::string cname);

    /**
     * This method search the materials map with the name string and return
     * a pointer to a Color object if it exists, otherwise it returns 0
     */
    Color* getMaterialColorByName(std::string cname);

    /**
     * This method search the numerical constant values map with the
     * name string and return the value if it exists,otherwhise it
     * returns 0 */
    double getNumConstByName(std::string cname);
    /**
     * This method search the char constant values map with the name
     * string and return the value if it exists,otherwhise it returns
     * 0 */
    std::string getCharConstByName(std::string cname);

    /** This method return the total number of sections contained 
	in the XML file */
    unsigned int getSectionsNumber();
    /** This method return the total number of volumes in the XML file */
    unsigned int getVolumesNumber();
    /** This method return the total number of constants in the XML file */
    unsigned int getConstantsNumber();
    /** This method return the total number of constants in the XML file */
    unsigned int getMaterialsNumber();

    /// This methods return the ith constant in the constants map
    Const* getOrderedConst(int i);
    /// This methods return the ith volume in the volumes map
    Volume* getOrderedVolume(int i);

    /// Method to retrive the finest mode string from the section
    std::string getFineChoice(){return (sections[0])->getFineChoice();}

    /** This method resolve the symbolic volumes references in the XML file 
     *  and it is automatically called by the builder.
     */
    void ResolveReferences();

  private:
    /** @link aggregation */
    /// This is the sections container
    std::vector < Section * > sections;  
  
    /** This is a private pointer to the Gddconstants that provide 
	the entry point to all the generic model hierarchy of constants*/
    Constants * constants;

    /// This is the collection of materials
    MatCollection* materials;

    /// Representation of Id dictionary 
    IdDict* dict;

    /** This is a private Gddvolume map indicized by names */
    std::map < std::string, Volume * > volumeMap;
    /** This is a private Gddconst map indicized by names */
    std::map < std::string, Const * > constMap;
    /** This is a private map of numerical constants values indicized by names */
    std::map < std::string, double > constNumMap;
    /** This is a private map of string constants values indicized by names */
    std::map < std::string, std::string > constCharMap;
    /// This is a vector of strings containing the possible modes
    std::vector <std::string> modeNames;

    /// This is the DTDversion
    std::string DTDversion;
    /// This is the CVSid
    std::string CVSid;
  };
}

#endif //GDD_H







#ifndef GDD_H
#define GDD_H
#include <map.h>
#include <vector.h>
#include <string.h>
#include <list>

class GDDvolume;
class GDDconstants;
class GDDconst;
class GDDvisitor;
class GDDsection;
class GDDchoice;
class GDDshape;

typedef vector <GDDsection*> vec;

/**
 * This is the main container of all the geometry.
 * Responsability of its creation and destruction is of the manager.
 * This class provides a series of maps and methods to query that 
 * maps. 
 * 
 * @author D.Favretto & R.Giannitrapani
 */
class GDD {
public:
  /** This is the destructor; it should be called only by the manager
   * destructor. It starts the destruction of all the objects created
   * by the builder
   */
  virtual ~GDD();

  /// This method gives back the DTDversion
  string getDTDversion(){return DTDversion;};
  /// This method gives back the CVSid  
  string getCVSid(){return CVSid;};
  /// This method sets the DTDversion  
  void setDTDversion(string pdtd){DTDversion = pdtd;};
  /// This method sets the CVSid  
  void setCVSid(string pcvs){CVSid = pcvs;};

  ///This method sets the constants 
  void setConstants(GDDconstants* pconstants){constants=pconstants;}
  ///This method gives back the constants
  GDDconstants* getConstants()const{return constants;}
  /// This method gives back a pointer to the sections vector
  vec* getSections(){return &sections;};  
  ///This method gives back the GDDvolume map  
  map < string, GDDvolume * > getVolumesMap(){return volumeMap;};

  /// This is the recursive accept for the visitor pattern
  void Accept(GDDvisitor* v);
  /// This is the non recursive accept for the visitor pattern
  void AcceptNotRec(GDDvisitor* v);

  /// This method build a global choices map for all the sections 
  void buildChoiceMap();
  /// This method build a global volumes map for all the sections 
  void buildVolumeMap();
  /// This method build the global constants maps for all the constants
  void buildConstantsMap();

  void buildBoundingBoxes(); 
  /// This method gives back the modes names vector
  vector <string> getModeNames(){return modeNames;};
  /// This method gives back the materials names vector
  vector <string> getMaterialNames(){return materialNames;};

  /**
   * This method search the volumes map with the name string and return 
   * a pointer to the GDDvolume if it exists, otherwise it returns a null 
   * pointer. It uses also the GetChoiceByName to search for possible choices
   * depending on the mode (see the manager).
   */
  GDDvolume * getVolumeByName(string vname);

  GDDshape * getShapeByName(string vname);
  /**
   * This method search the const map with the name string and return
   * a pointer to the GDDconst object if it exists, otherwise it returns 0
   */
  GDDconst* getConstByName(string cname);
  /**
   * This method search the numerical constant values map with the
   * name string and return the value if it exists,otherwhise it
   * returns 0 */
  double getNumConstByName(string cname);
  /**
   * This method search the char constant values map with the name
   * string and return the value if it exists,otherwhise it returns
   * 0 */
  string getCharConstByName(string cname);
  /**
   * This method search the choices map with the name string and return 
   * a pointer to the GDDchoice if it exists, otherwise it returns a null 
   * pointer. 
   */  
  GDDchoice * getChoiceByName(string cname);

  /** This method return the total number of sections contained 
      in the XML file */
  int getSectionsNumber();
  /** This method return the total number of volumes in the XML file */
  int getVolumesNumber();


  /** This method resolve the symbolic volumes references in the XML file 
   *  and it is automatically called by the builder.
   */
  void ResolveReferences();

  
  void visitMap();

 private:
  /** @link aggregation */
  /// This is the sections container
  vector < GDDsection * > sections;  
  
  /** This is a private pointer to the GDDconstants that provide 
      the entry point to all the generic model hierarchy of constants*/
  GDDconstants * constants;

  /** This is a private GDDchoice map indicized by names */
  map < string, GDDchoice * > choiceMap;  
  /** This is a private GDDvolume map indicized by names */
  map < string, GDDvolume * > volumeMap;
  /** This is a private GDDconst map indicized by names */
  map < string, GDDconst * > constMap;
  /** This is a private map of numerical constants values indicized by names */
  map < string, double > constNumMap;
  /** This is a private map of string constants values indicized by names */
  map < string, string > constCharMap;

  /// This is a vector of strings containing the possible modes
  vector <string> modeNames;
  /// This is a vector of strings containing the material names
  vector <string> materialNames;

  /// This is the DTDversion
  string DTDversion;
  /// This is the CVSid
  string CVSid;
};
#endif //GDD_H







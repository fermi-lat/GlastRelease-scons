#ifndef GDDVOLUME_H
#define GDDVOLUME_H
#include <string>

class GDDsectionsVisitor;

/// All the enumeratives that are needed in the volumes hierarchy
enum unitLengthtype{GDDmm,GDDcm,GDDm};
enum unitAngletype{GDDdeg,GDDmrad};
enum stacktype{sx,sy,sz};
enum stackorigin{atStart,atCentre};
enum typeVolume{shape,composition,stack,logical};

/**
 * This class represent the generic volume of the geometric 
 * description. It is not an abstract class, but it should be not
 * explicitely constructed. Its purpouse is to simplify the volumes
 * hierarchy in the GDD.  
 * @author D.Favretto & R.Giannitrapani */
class GDDvolume {
 public:
  /**
   * This method sets the parameters associated with the volume
   */
  void setParameters(string pparameters){parameters = pparameters;};
  /**
   * This method sets the volume name 
   */
  void setName(string pname){name = pname;};
  /**
   * This method returns a string with the name of the volume
   */
  string getName(){return name;};
  /**
   * This method returns a string with the parameters
   */
  string getParameters(){return parameters;};
  /**
   * This method returns the type of the volume
   */  
  typeVolume getVolumeType(){return vtype;};

  /** This method return the x dimension of the bounding box it has to
     be redefined concretely in the subclasses of GDDvolume */
  virtual double getBBX(){return 0;};
  /** This method return the y dimension of the bounding box it has to
     be redefined concretely in the subclasses of GDDvolume */
  virtual double getBBY(){return 0;};
  /** This method return the z dimension of the bounding box it has to
     be redefined concretely in the subclasses of GDDvolume */
  virtual double getBBZ(){return 0;};
  
  /// This method is the recursive Accept for the visitor pattern
  virtual void Accept(GDDsectionsVisitor*){};
  /// This method is the non recursive Accept for the visitor pattern
  virtual void AcceptNotRec(GDDsectionsVisitor*){};
  
  virtual ~GDDvolume(){}

 protected:

  /// Some constructors
  GDDvolume(){};
  GDDvolume(string pname,typeVolume ptype):
    name(pname),vtype(ptype),parameters(""){;}
  GDDvolume(typeVolume ptype):vtype(ptype),parameters(""){;}



 private:
  /// The name of the volume
  string name;
  /// The type of the volumes
  typeVolume vtype;
  /// The parameters of the volume
  string parameters;  

};
#endif //GDDVOLUME_H




#ifndef GDDVOLUME_H
#define GDDVOLUME_H
#include <string>

class GDDsectionsVisitor;

/// All the enumeratives that are needed in the volumes hierarchy
enum unitLengthtype{GDDmm,GDDcm,GDDm};
enum unitAngletype{GDDdeg,GDDmrad};
enum typeVolume{shape,composition,stack,logical};
enum stacktype{sx,sy,sz};


/**
 * This class represent the generic volume of the geometric 
 * description. It is not an abstract class, but it should be not
 * explicitely constructed. Its purpouse is to simplify the volumes
 * hierarchy in the GDD.  
 * @author D.Favretto & R.Giannitrapani */
class GDDvolume {
 public:
  /// Some constructors
  GDDvolume(){};
  GDDvolume(std::string pname,typeVolume ptype):
    name(pname),vtype(ptype),parameters(""),bbx(0),bby(0),bbz(0){;}
  GDDvolume(typeVolume ptype):vtype(ptype),parameters(""){;}

  /**
   * This method sets the parameters associated with the volume
   */
  void setParameters(std::string pparameters){parameters = pparameters;};
  /**
   * This method sets the volume name 
   */
  void setName(std::string pname){name = pname;};
  /**
   * This method returns a string with the name of the volume
   */
  std::string getName(){return name;};
  /**
   * This method returns a string with the parameters
   */
  std::string getParameters(){return parameters;};
  /**
   * This method returns the type of the volume
   */  
  typeVolume getVolumeType(){return vtype;};

  
  /// This method is the recursive Accept for the visitor pattern
  virtual void Accept(GDDsectionsVisitor*) = 0;
  /// This method is the non recursive Accept for the visitor pattern
  virtual void AcceptNotRec(GDDsectionsVisitor*) = 0;
  
  /**
   * This method build the bounding box of the volume.
   * It is concretely implemented in subclasses of GDDvolume
   */
  virtual void constructBB() = 0;

  /** This method return the x dimension of the bounding box it has to
      be redefined concretely in the subclasses of GDDvolume */
  double getBBX(){return bbx;};
  /** This method return the y dimension of the bounding box it has to
     be redefined concretely in the subclasses of GDDvolume */
  double getBBY(){return bby;};
  /** This method return the z dimension of the bounding box it has to
     be redefined concretely in the subclasses of GDDvolume */
  double getBBZ(){return bbz;};

  /// This method sets the x dimension of the bounding box 
  void setBBX(double pbbx){bbx = pbbx;};
  /// This method sets the y dimension of the bounding box 
  void setBBY(double pbby){bby = pbby;};
  /// This method sets the z dimension of the bounding box 
  void setBBZ(double pbbz){bbz = pbbz;};

  /** This method return the x dimension of the bounding box it has to
      be redefined concretely in the subclasses of GDDvolume */
  double getOX(){return ox;};
  /** This method return the y dimension of the bounding box it has to
     be redefined concretely in the subclasses of GDDvolume */
  double getOY(){return oy;};
  /** This method return the z dimension of the bounding box it has to
     be redefined concretely in the subclasses of GDDvolume */
  double getOZ(){return oz;};

  /// This method sets the x dimension of the bounding box 
  void setOX(double pox){ox = pox;};
  /// This method sets the y dimension of the bounding box 
  void setOY(double poy){oy = poy;};
  /// This method sets the z dimension of the bounding box 
  void setOZ(double poz){oz = poz;};



  
  virtual ~GDDvolume(){}

 private:

  /// The name of the volume
  std::string name;
  /// The type of the volumes
  typeVolume vtype;
  /// The parameters of the volume
  std::string parameters;  
  /// Bounding box
  double bbx;
  double bby;
  double bbz;

  double ox;
  double oy;
  double oz;
};
#endif //GDDVOLUME_H




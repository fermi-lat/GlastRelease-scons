#ifndef GDDVOLUME_H
#define GDDVOLUME_H

#include "detModel/Sections/GDDboundable.h"
#include <string>

class GDDsectionsVisitor;

/**
 * This class represent the generic volume of the geometric 
 * description. It is an abstract class. Since a volume has 
 * a bounding box, this class is a subclass of GDDboundable
 * @author R.Giannitrapani & D.Favretto
 */
class GDDvolume: public GDDboundable {

 public:
  /** A constructors with type and name */
  GDDvolume(std::string pname):name(pname),parameters(""){;}
  /** A constructors with only type */
  GDDvolume():name(""),parameters(""){;}

  ~GDDvolume(){;};

  /**
   * This method sets and gets the parameters associated with the
   * volume.  Parameters can be any string, representing a way to
   * extend the class in the future.  
   */
  void setParameters(std::string pparameters){parameters = pparameters;};
  std::string getParameters(){return parameters;};
  /**
   * This methods sets and gets the volume name 
   */
  void setName(std::string pname){name = pname;};
  std::string getName(){return name;};

  /// This method is the recursive Accept for the visitor pattern
  virtual void Accept(GDDsectionsVisitor*) = 0;
  /// This method is the non recursive Accept for the visitor pattern
  virtual void AcceptNotRec(GDDsectionsVisitor*) = 0;

 private:

  /// The name of the volume
  std::string name;
  /// The parameters of the volume
  std::string parameters;  
};
#endif //GDDVOLUME_H




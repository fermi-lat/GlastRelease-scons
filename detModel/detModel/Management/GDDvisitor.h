#ifndef GDDVISITOR_H
#define GDDVISITOR_H
#include <string>

/// This enumerative represent the type of visitor
enum typeVisitor{sectionsVisitor, constantsVisitor};

/**
 * This abstract class is the base for the visitors hierarchy; its
 * main purpouse is to provide an abstract interface common to all the
 * different kind of visitors (one type for each packages); in this
 * way the manager can hold a list of visitors that the client can
 * start.  
 * @author D.Favretto and R.Giannitrapani 
*/
class GDDvisitor {

 public:

  /// This method gives back the type of the visitor
  typeVisitor getType(){return type;};
  /// This method sets the type of the visitor
  void setType(typeVisitor ptype){type = ptype;};
  /// This method says if the visitor is recursive or not
  bool getRecursive(){return recursive;}
  /// This method sets if the visitor is recursive or not
  void setRecursive(bool prec){recursive = prec;};

 protected:
  GDDvisitor(){};

  /// If this is true, the visitor is recursive, otherwise it is not.
  bool recursive;
  /**
   * This is the type of the visitor; possible types are listed in the
   * enumerative type typeVisitor  */
  typeVisitor type;
};
#endif //GDDVISITOR_H

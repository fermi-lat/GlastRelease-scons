#ifndef VISITOR_H
#define VISITOR_H
#include <string>

namespace detModel{

  class Gdd;

  /**
   * This abstract class is the base for the visitors hierarchy; its
   * main purpouse is to provide an abstract interface common to all the
   * different kind of visitors (one type for each packages); in this
   * way the manager can hold a list of visitors that the client can
   * start.  
   * @author D.Favretto and R.Giannitrapani 
   */
  class Visitor {
  public:
    /// This method says if the visitor is recursive or not
    bool getRecursive(){return recursive;}
    /// This method sets if the visitor is recursive or not
    void setRecursive(bool prec){recursive = prec;};
    
    virtual void visitGdd(Gdd*) = 0;
    
  protected:
    /// If this is true, the visitor is recursive, otherwise it is not.
    bool recursive;
  };
}
#endif //VISITOR_H

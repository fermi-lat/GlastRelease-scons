#ifndef BUILDER_H
#define BUILDER_H
#include "detModel/Management/Manager.h"

namespace detModel{


  class Gdd;
  class Constants;


  /**
   * This class defines the public interface for the generic
   * functionalities of a builder. Each subclass must provide
   * implementation for some of its methods. For this it is not a pure
   * abstract class.  @author D.Favretto & R.Giannitrapani 
   */
  class Builder {
  public:
  
    Builder(){};

    /** This method initialize the parser (implemented in a concrete
     * subclass of this class) 
     */
    virtual void parseFile(char*) = 0;

    /** This method start the parser (implemented in a concrete
     * subclass of this class) for the sections part and put the
     * result in a private pointer */
    virtual void buildSections() = 0;

    /**
     * This method start the parser (implemented in a concrete
     * subclass of this class) for the constants part and put the
     * result in a private pointer */
    virtual void buildConstants() = 0;

    void setGdd(Gdd* pGdd){currentGdd = pGdd;}

    /** This method return a pointer to the internal
        Gdd                  
    */
    Gdd * getGdd() { return currentGdd; };

  protected:
    /** This is a private pointer to a Gdd used
        during the building */
    Gdd * currentGdd;
  };

}
#endif //BUILDER_H


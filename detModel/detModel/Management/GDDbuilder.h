#ifndef GDDBUILDER_H
#define GDDBUILDER_H
#include "detModel/Management/GDDmanager.h"

class GDD;
class GDDconstants;


/**
 * This class defines the public interface for the generic functionalities
 * of a builder. Each subclass must provide implementation for some of its methods. For this it is not a pure abstract class.
 * @author D.Favretto & R.Giannitrapani
 */
class GDDbuilder {
public:

    /** This method start the parser (implemented in a concrete
     * subclass of this class) for the sections part and put the
     * result in a private pointer */
    virtual void buildSections() { };

    /**
     * This method start the parser (implemented in a concrete
     * subclass of this class) for the constants part and put the
     * result in a private pointer */
    virtual void buildConstants() { };

    /** This method return a pointer to the internal
        GDDsectionsContainer */
    GDD * getGDD() { return currentGDD; };

protected:
    GDDbuilder(){};

    /** This is a private pointer to a GDDsectionsContainer used
        during the building */
    GDD * currentGDD;
};
#endif //GDDBUILDER_H

#ifndef VRMLSECTIONSVISITOR_H
#define VRMLSECTIONSVISITOR_H
#include "detModel/Management/SectionsVisitor.h"
#include <fstream>
#include <vector>
#include <map>

namespace detModel{

  class Ensemble;
  class Gdd;
  class Color;

  /*
   * This is a concrete implementation of a sectionsVisitor that produces
   * a VRML file with the geometry. 
   */
  class VrmlSectionsVisitor : public SectionsVisitor {

  public:

    /**
     * The constructor has a string parameter; it is the name of the
     * root volume in GDD to represent along with all the volumes
     * contained in it. If it is the empty string the root volume is the
     * topVolume of the section.  */
    VrmlSectionsVisitor(std::string); 
  
    /// A standard destructor
    virtual ~VrmlSectionsVisitor();
  
    /**
     * This is the visitor for the GDDsectionsContainer */
    virtual void visitGdd(Gdd*);
  
    /**
     * This is the visitor for the Section 
     */
    virtual void visitSection(Section*);

    /**
     * This is the visitor for the Ensemble 
     */
    virtual void visitEnsemble(Ensemble*);

    /**
     * This is the visitor for the Box 
     */
    virtual void visitBox(Box*);

    /**
     * This is the visitor for the Box 
     */
    virtual void visitTube(Tube*);

    /**
     * This is the visitor for the PosXYZ 
     */
    virtual void visitPosXYZ(PosXYZ*);

    /**
     * This is the visitor for the AxisMPos 
     */
    virtual void visitAxisMPos(AxisMPos*);

    /**
     * This is the visitor for the IdField 
     */
    virtual void visitIdField(IdField*);

    /**
     * This is the visitor for the Seg 
     */
    virtual void visitSeg(Seg*);

    /// This method build the colors for the VRML file
    void makeColor();

    void setDepth(std::string, int);
  private:
    std::string actualVolume;

    std::ofstream out;

    /** This variables keep the actual depth or recursivness of the
	visit in the hierarchy */
    int depth;
    /// This map holds the colors for the material
    std::map <std::string, Color*> colorsMap;
    /// This map holds information on the recursion depth of volumes
    std::map <std::string, int> depthMap;

    int numVol;
    Gdd* m_gdd;
  };

}
#endif //VRMLSECTIONSVISITOR_H








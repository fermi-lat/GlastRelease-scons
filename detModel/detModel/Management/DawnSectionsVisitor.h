#ifndef DAWNSECTIONSVISITOR_H
#define DAWNSECTIONSVISITOR_H
#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Utilities/Vector.h"
#include <fstream>
#include <vector>
#include <map>


namespace detModel{

  class Gdd;
  class Color;
  /*
   * This is a concrete implementation of a sectionsVisitor that produces
   * a DAWN file with the geometry. 
   */
  class DawnSectionsVisitor : public SectionsVisitor {
    
  public:
    
    DawnSectionsVisitor(std::string);
    virtual ~DawnSectionsVisitor();
    
    /**
     * This is the visitor for the SectionsContainer 
     */
    virtual void visitGdd(Gdd*);
    
    /**
     * This is the visitor for the Section 
     */
    virtual void visitSection(Section*);
    
    /**
     * This is the visitor for the GDDcomposition 
     */
    virtual void visitEnsemble(Ensemble*);
    
    /**
     * This is the visitor for the GDDbox 
     */
    virtual void visitBox(Box*);

    /**
     * This is the visitor for the Tube 
     */
    virtual void visitTube(Tube*){};
    
    /**
     * This is the visitor for the GDDposXYZ 
     */
    virtual void visitPosXYZ(PosXYZ*);
    
    /**
     * This is the visitor for the GDDaxisMPos 
     */
    virtual void visitAxisMPos(AxisMPos*);
    
    /**
     * This is the visitor for the GDDaxisMPos 
     */
    virtual void visitIdField(IdField*);
    
    /**
     * This is the visitor for the GDDseg 
     */
    virtual void visitSeg(Seg*);
    
    std::string actualVolume;
    
    std::ofstream out;
    
    Vector* actualPos;
    Vector* actualXAxis;
    Vector* actualYAxis;
    Vector* actualZAxis;
    
    double actualXRot;
    double actualYRot;
    double actualZRot;

    /// This map holds the colors for the material
    std::map <std::string, Color*> colorsMap;
    
  };
}
#endif //DAWNSECTIONSVISITOR_H








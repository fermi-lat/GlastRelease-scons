#ifndef GDDVRMLSECTIONSVISITOR_H
#define GDDVRMLSECTIONSVISITOR_H
#include "detModel/Management/GDDsectionsVisitor.h"
#include <fstream>
#include <vector>
#include <map>

class GDD;

/*
 * This is a concrete implementation of a sectionsVisitor that produces
 * a VRML file with the geometry. 
 */
class GDDVRMLSectionsVisitor : public GDDsectionsVisitor {

 public:

  /**
   * The constructor has a string parameter; it is the name of the
   * root volume in GDD to represent along with all the volumes
   * contained in it. If it is the empty string the root volume is the
   * topVolume of the section.  */
  GDDVRMLSectionsVisitor(std::string); 
  
  /// A standard destructor
  virtual ~GDDVRMLSectionsVisitor();
  
  /**
   * This is the visitor for the GDDsectionsContainer */
  virtual void visitGDD(GDD*);
  
  /**
   * This is the visitor for the GDDsection 
   */
  virtual void visitSection(GDDsection*);

  /**
   * This is the visitor for the GDDensamble 
   */
  virtual void visitEnsamble(GDDensamble*);

  /**
   * This is the visitor for the GDDbox 
   */
  virtual void visitBox(GDDbox*);

  /**
   * This is the visitor for the GDDposXYZ 
   */
  virtual void visitPosXYZ(GDDposXYZ*);

  /**
   * This is the visitor for the GDDaxisPos 
   */
  virtual void visitAxisPos(GDDaxisPos*);

  /**
   * This is the visitor for the GDDaxisMPos 
   */
  virtual void visitAxisMPos(GDDaxisMPos*);

  /**
   * This is the visitor for the GDDaxisMPos 
   */
  virtual void visitIdField(GDDidField*);

  /**
   * This is the visitor for the GDDseg 
   */
  virtual void visitSeg(GDDseg*);

  /// This method build the colors for the VRML file
  void makeColor();

  void setOpacity(std::string, float);

  void setDepth(std::string, int);

  std::string actualVolume;

  std::ofstream out;

  /** This variables keep the actual depth or recursivness of the
      visit in the hierarchy */
  int depth;
  /// This map holds the opacity information of the material colors
  std::map <std::string, float> opacityMap;
  /// This map holds information on the recursion depth of volumes
  std::map <std::string, int> depthMap;
};
#endif //GDDVRMLSECTIONSVISITOR_H








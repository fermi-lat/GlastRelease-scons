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

  GDDVRMLSectionsVisitor(string);
  virtual ~GDDVRMLSectionsVisitor();
  
  /**
   * This is the visitor for the GDDsectionsContainer 
   */
  virtual void visitGDD(GDD*);
  
  /**
   * This is the visitor for the GDDsection 
   */
  virtual void visitSection(GDDsection*);

  /**
   * This is the visitor for the GDDcomposition 
   */
  virtual void visitComposition(GDDcomposition*);

  /**
   * This is the visitor for the GDDstack 
   */
  virtual void visitStack(GDDstack*);

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

  string actualVolume;

  ofstream out;

  /// This map holds the opacity information of the material colors
  map <string, float> opacityMap;
  /// This map holds information on the recursion depth of volumes
  map <string, int> depthMap;
};
#endif //GDDVRMLSECTIONSVISITOR_H








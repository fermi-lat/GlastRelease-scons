#ifndef GDDVRMLSECTIONSVISITOR_H
#define GDDVRMLSECTIONSVISITOR_H
#include "detModel/Management/GDDsectionsVisitor.h"

class GDD;

/*
 * This is a concrete implementation of a sectionsVisitor that produces
 * a VRML file with the geometry. 
 */
class GDDVRMLSectionsVisitor : public GDDsectionsVisitor {

 public:

  GDDVRMLSectionsVisitor();
  GDDVRMLSectionsVisitor(string);
  
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

};
#endif //GDDVRMLSECTIONSVISITOR_H








#ifndef GDDPRINTERSECTIONSVISITOR_H
#define GDDPRINTERSECTIONSVISITOR_H
#include "detModel/Management/GDDsectionsVisitor.h"
class GDD;

/*
 * This is a concrete implementation of a sectionsVisitor. It
 * simply traverses all the sections hierarchy and print on the
 * standard output some information on the geometry. It is
 * a recursive visitor.
 */
class GDDprinterSectionsVisitor : public GDDsectionsVisitor {

 public:
  /// This is the constructor
  GDDprinterSectionsVisitor();
  
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
};
#endif //GDDPRINTERSECTIONSVISITOR_H








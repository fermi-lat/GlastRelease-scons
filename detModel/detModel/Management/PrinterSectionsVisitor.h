#ifndef PRINTERSECTIONSVISITOR_H
#define PRINTERSECTIONSVISITOR_H
#include "detModel/Management/SectionsVisitor.h"

namespace detModel{

  /*
   * This is a concrete implementation of a sectionsVisitor. It
   * simply traverses all the sections hierarchy and print on the
   * standard output some information on the geometry. It is
   * a recursive visitor.
   */
  class PrinterSectionsVisitor : public SectionsVisitor {

  public:
    /// This is the constructor
    PrinterSectionsVisitor();
  
    /**
     * This is the visitor for the Gdd 
     */
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
  };

}
#endif //PRINTERSECTIONSVISITOR_H








#ifndef GDDPRINTERCONSTANTSVISITOR_H
#define GDDPRINTERCONSTANTSVISITOR_H
#include "detModel/Management/GDDconstantsVisitor.h"

/*
 * This is a concrete implementation of a constantsVisitor. It
 * simply traverses all the constants hierarchy and print on the
 * standard output some information on them. It is
 * a recursive visitor.
 */
class GDDprinterConstantsVisitor : public GDDconstantsVisitor {

 public:
  /// This is the constructor
  GDDprinterConstantsVisitor();
  
  /**
   * This is the visitor for the GDD
   */
  virtual void visitGDD(GDD*);
  /**
   * This is the visitor for the GDDconstants
   */
  virtual void visitConstants(GDDconstants*);
  /**
   * This is the visitor for the GDDconstCategory
   */
  virtual void visitConstCategory(GDDconstCategory*);
  /**
   * This is the visitor for the GDDintConst 
   */
  virtual void visitIntConst(GDDintConst*);
  /**
   * This is the visitor for the GDDfloatConst 
   */
  virtual void visitFloatConst(GDDfloatConst*);
  /**
   * This is the visitor for the GDDdoubleConst 
   */
  virtual void visitDoubleConst(GDDdoubleConst*);
  /**
   * This is the visitor for the GDDstringConst 
   */
  virtual void visitStringConst(GDDstringConst*);
};
#endif // GDDPRINTERCONSTANTSVISITOR_H









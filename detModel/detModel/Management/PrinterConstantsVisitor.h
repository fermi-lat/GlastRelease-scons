#ifndef PRINTERCONSTANTSVISITOR_H
#define PRINTERCONSTANTSVISITOR_H
#include "detModel/Management/ConstantsVisitor.h"

namespace detModel{

  /*
   * This is a concrete implementation of a constantsVisitor. It
   * simply traverses all the constants hierarchy and print on the
   * standard output some information on them. It is
   * a recursive visitor.
   */
  class PrinterConstantsVisitor : public ConstantsVisitor {

  public:
    /// This is the constructor
    PrinterConstantsVisitor();
  
    /**
     * This is the visitor for the Gdd
     */
    virtual void visitGdd(Gdd*);
    /**
     * This is the visitor for the Constants
     */
    virtual void visitConstants(Constants*);
    /**
     * This is the visitor for the ConstCategory
     */
    virtual void visitConstCategory(ConstCategory*);
    /**
     * This is the visitor for the IntConst 
     */
    virtual void visitIntConst(IntConst*);
    /**
     * This is the visitor for the FloatConst 
     */
    virtual void visitFloatConst(FloatConst*);
    /**
     * This is the visitor for the DoubleConst 
     */
    virtual void visitDoubleConst(DoubleConst*);
    /**
     * This is the visitor for the StringConst 
     */
    virtual void visitStringConst(StringConst*);
  };

}
#endif // PRINTERCONSTANTSVISITOR_H









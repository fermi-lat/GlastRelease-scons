#ifndef CONSTANTSVISITOR_H
#define CONSTANTSVISITOR_H
#include "detModel/Management/Visitor.h"

namespace detModel{

  class Gdd;
  class Constants;
  class ConstCategory;
  class IntConst;
  class FloatConst;
  class DoubleConst;
  class StringConst;

  /**
   * This is the abstract class for constants visitors, i.e. visitors
   * that have to traverse the constants of the XML file */
  class ConstantsVisitor : public Visitor {
  public:
    /**
     * This is the visitor for the Gdd
     */
    virtual void visitGdd(Gdd*) = 0;
    /**
     * This is the visitor for the Constants
     */
    virtual void visitConstants(Constants*) = 0;
    /**
     * This is the visitor for the ConstCategory
     */
    virtual void visitConstCategory(ConstCategory*) = 0;
    /**
     * This is the visitor for the IntConst 
     */
    virtual void visitIntConst(IntConst*) = 0;
    /**
     * This is the visitor for the FloatConst 
     */
    virtual void visitFloatConst(FloatConst*) = 0;
    /**
     * This is the visitor for the DoubleConst 
     */
    virtual void visitDoubleConst(DoubleConst*) = 0;
    /**
     * This is the visitor for the StringConst 
     */
    virtual void visitStringConst(StringConst*) = 0;
  };
}

#endif //CONSTANTSVISITOR_H


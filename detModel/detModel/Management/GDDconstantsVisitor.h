#ifndef GDDCONSTANTSVISITOR_H
#define GDDCONSTANTSVISITOR_H
#include "detModel/Management/GDDvisitor.h"

class GDD;
class GDDconstants;
class GDDconstCategory;
class GDDintConst;
class GDDfloatConst;
class GDDdoubleConst;
class GDDstringConst;

/**
 * This is the abstract class for constants visitors, i.e. visitors that have to traverse the constants of the XML file 
 */
class GDDconstantsVisitor : public GDDvisitor {
 public:
  /**
   * This is the visitor for the GDD
   */
  virtual void visitGDD(GDD*) = 0;
  /**
   * This is the visitor for the GDDconstants
   */
  virtual void visitConstants(GDDconstants*) = 0;
  /**
   * This is the visitor for the GDDconstCategory
   */
  virtual void visitConstCategory(GDDconstCategory*) = 0;
  /**
   * This is the visitor for the GDDintConst 
   */
  virtual void visitIntConst(GDDintConst*) = 0;
  /**
   * This is the visitor for the GDDfloatConst 
   */
  virtual void visitFloatConst(GDDfloatConst*) = 0;
  /**
   * This is the visitor for the GDDdoubleConst 
   */
  virtual void visitDoubleConst(GDDdoubleConst*) = 0;
  /**
   * This is the visitor for the GDDstringConst 
   */
  virtual void visitStringConst(GDDstringConst*) = 0;
};
#endif //GDDCONSTANTSVISITOR_H


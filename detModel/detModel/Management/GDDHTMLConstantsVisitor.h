#ifndef GDDHTMLCONSTANTSVISITOR_H
#define GDDHTMLCONSTANTSVISITOR_H
#include "detModel/Management/GDDconstantsVisitor.h"
#include <fstream>

/*
 * This is a concrete implementation of a constantsVisitor. It
 * simply traverses all the constants hierarchy and print on a
 * html output some information on them. It is
 * a recursive visitor.
 */
class GDDHTMLConstantsVisitor : public GDDconstantsVisitor {

 public:
  /// This is the constructor
  GDDHTMLConstantsVisitor();
  /// This is the destructor
  virtual ~GDDHTMLConstantsVisitor();
  
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
   * This is the visitor for the GDDconst 
   */
  virtual void visitConst(GDDconst*);

 private:
  std::ofstream out;

};
#endif // GDDHTMLCONSTANTSVISITOR_H









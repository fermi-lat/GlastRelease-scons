#ifndef PRINTERMATERIALSVISITOR_H
#define PRINTERMATERIALSVISITOR_H
#include "detModel/Management/MaterialsVisitor.h"

namespace detModel{

  class Gdd;
  class MatCollection;
  class Element;
  class Composite;

  /**
   * This is the abstract class for materials visitors, i.e. visitors
   * that have to traverse the materials of the XML file 
   * @author M.Frailis & R.Giannitrapani 
   */

  class PrinterMaterialsVisitor : public MaterialsVisitor {
  public:
    PrinterMaterialsVisitor();

    /**
     * This is the visitor for the Gdd
     */
    virtual void visitGdd(Gdd*);
    /**
     * This is the visitor for the MatCollection
     */
    virtual void visitMatCollection(MatCollection*);
    /**
     * This is the visitor for the Element
     */
    virtual void visitElement(Element*);
    /**
     * This is the visitor for the Composite 
     */
    virtual void visitComposite(Composite*);
  };
}

#endif //PRINTERMATERIALSVISITOR_H


#ifndef MATERIALSVISITOR_H
#define MATERIALSVISITOR_H
#include "detModel/Management/Visitor.h"

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

  class MaterialsVisitor : public Visitor {
  public:
    /**
     * This is the visitor for the Gdd
     */
    virtual void visitGdd(Gdd*) = 0;
    /**
     * This is the visitor for the MatCollection
     */
    virtual void visitMatCollection(MatCollection*) = 0;
    /**
     * This is the visitor for the Element
     */
    virtual void visitElement(Element*) = 0;
    /**
     * This is the visitor for the Composite 
     */
    virtual void visitComposite(Composite*) = 0;
  };
}

#endif //MATERIALSVISITOR_H


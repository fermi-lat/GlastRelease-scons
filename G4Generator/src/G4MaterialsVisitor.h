#ifndef G4MATERIALSVISITOR_H
#define G4MATERIALSVISITOR_H
#include "detModel/Management/MaterialsVisitor.h"
#include <fstream>
#include <vector>
#include <map>

class G4MaterialsVisitor: public detModel::MaterialsVisitor {

 public:

  G4MaterialsVisitor();
  virtual ~G4MaterialsVisitor(){};
  
  /**
   * This is the visitor for the Gdd
   */
  virtual void visitGdd(detModel::Gdd*);
  /**
   * This is the visitor for the MatCollection
   */
  virtual void visitMatCollection(detModel::MatCollection*);
  /**
   * This is the visitor for the Element
   */
  virtual void visitElement(detModel::Element*);
  /**
   * This is the visitor for the Composite 
   */
  virtual void visitComposite(detModel::Composite*);
};

#endif //G4MATERIALSVISITOR_H






#ifndef MEDIAVISITOR_H
#define MEDIAVISITOR_H
#include "detModel/Management/MaterialsVisitor.h"
#include "GlastSvc/GlastDetSvc/IMedia.h"
#include <fstream>
#include <vector>

/** @class MediaVisitor
* @brief
*
* $Header$
*/
class MediaVisitor: public detModel::MaterialsVisitor {

 public:

  MediaVisitor(IMedia&);
  virtual ~MediaVisitor(){};
  
  /**
   * This is the visitor for the Gdd
   */
  virtual void visitGdd(detModel::Gdd*){};
  /**
   * This is the visitor for the MatCollection
   */
  virtual void visitMatCollection(detModel::MatCollection*){};
  /**
   * This is the visitor for the Element
   */
  virtual void visitElement(detModel::Element*);
  /**
   * This is the visitor for the Composite 
   */
  virtual void visitComposite(detModel::Composite*);

 private:
  /// The list of materials names
  std::vector <std::string> m_materials;  

  //! the object to call back, passed in with the constructor
  IMedia& m_media;
};

#endif //MEDIAVISITOR_H

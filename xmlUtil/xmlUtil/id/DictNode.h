// $Header$
#ifndef XMLUTIL_DICTNODE_H
#define XMLUTIL_DICTNODE_H

#include "dom/DOM_Element.hpp"
#include "xmlUtil/id/DictField.h"
#include "xmlUtil/id/DictConstraints.h"

#include <vector>


namespace xmlUtil {
  //! DictNode represents a node in an identifier dictionary.
  /*! It could be embedded in IdDict, but it's too big and complicated
      for that to be very attractive.  Instead, make virtually everything
      private and let IdDict be a friend.
  */
  class DictNode {
  public:
    //! Construct a node from its XML representation (probably
    //! recursively construct children as well)
    DictNode(DOM_Element elt);
    ~DictNode();  // similarly, need to delete all children recursively

  private:
    friend class IdDict;
    friend class IdConverter;   //maybe

    // More constructor, for when dictionary is built dynamically
    // rather than from xml file
    DictNode(DictField& field);
    DictNode(DictField& field, DictConstraints& parConstraints);
    DictNode(DictField& field, DictConstraints& parConstraints,
             DictConstraints& myConstraints);
    const DictField& getField() const {return *field;};
    const Nodes&     getChildren() const {return children;};
    const DictConstraints& getValueConstraints() const { 
      return *myConstraints;};
    const DictConstraints& getParentConstraints() const {
      return *parConstraints;};
    bool  allowed(const unsigned value);
    bool  allowedChild(std::string childField, unsigned childValue, 
                       unsigned value);
    bool  allowedChild(unsigned childValue, unsigned value);
    bool  addChild(DictNode* child);

    typedef std::vector<DictNode*> Nodes;

    DictNode() {};                    /*< don't allow uninitialized node */
    Nodes            m_children;      /*< collection of child nodes */
    DictField       *m_field;         /*< id field for this node */
    DictNode        *m_parent;        /*< parent dictnode */
    DictConstraints *m_parConstraints;/*< Constraints on values of parent */
    DictConstraints *m_myConstraints; /*< Constraints on values for this node*/

  };       // end DictNode class definition
}          // end namespace xmlUtil
#endif

// $Header$
#ifndef XMLUTIL_DICTNODE_H
#define XMLUTIL_DICTNODE_H

#include "dom/DOM_Element.hpp"
#include "xmlUtil/id/DictField.h"
#include "xmlUtil/id/DictConstraints.h"

#include <vector>
#include <algorithm>
#include <set>
#include <iterator>

namespace xmlUtil {
  //! DictNode represents a node in an identifier dictionary.
  /*! It could be embedded in IdDict, but it's too big and complicated
      for that to be very attractive.  Instead, make virtually everything
      private and let IdDict be a friend.
  */
  class DictFieldMan;

  class DictNode {
  public:
    //! Construct a node from its XML representation (probably
    //! recursively construct children as well)
    DictNode(DOM_Element elt, DictNode* parent, DictFieldMan *fieldMan); 
    /* IN PROGRESS */
    ~DictNode();  // similarly, delete all children recursively /* TO DO */
    const DictField& getField() const {return *m_field;};
    const DictConstraints& getValueConstraints() const { 
      return *m_myConstraints;};
    const DictConstraints& getParentConstraints() const {
      return *m_parConstraints;};

  private:
    friend class IdDict;
    friend class IdConverter;   //maybe
    friend struct POrder;

    // Might want to use some other form of collection for 
    // child DictNodes
    typedef std::vector<DictNode*> Nodes;
    typedef std::vector<DictNode*>::iterator NodeIterator;

    // Nested class allowing sorting of child nodes by m_minValue of
    // parent constraint
    struct POrder {
      bool operator()(DictNode* a, DictNode* b) {
        return ((a->getParentConstraints()).getMin() <
                (b->getParentConstraints()).getMin()   ); }
    }; 

    const Nodes&     getChildren() const {return m_children;}
    // More constructor, for when dictionary is built dynamically
    // rather than from xml file
    /*    DictNode(DictField& field);                             TO DO
    DictNode(DictField& field, DictConstraints& parConstraints);  TO DO
    DictNode(DictField& field, DictConstraints& parConstraints,   TO DO
             DictConstraints& myConstraints);                     TO DO */
    bool  allowed(const unsigned value) const;            
    bool  allowedChild(std::string childField, unsigned childValue, 
                       unsigned value) const {return true;}         /* TO DO */
    bool  allowedChild(unsigned childValue, unsigned value) const {
      return true;  }  /* TO DO */
    bool  addChild(DictNode* child)    {return true;}            /* TO DO */


    //! Check that collection of children define appropriately disjoint
    //! collections of Identifiers.  In the process may do some sorting
    //! of child nodes, so these functions are not const.  
    bool  consistent();                        /* done, sort of */

    //! Insert allowed values for this DictNode into a set so we
    //! can use set operations to look for overlap
    void insertValues(std::set<unsigned>& aSet);

    //! Check that for nodes in the interval demarcated by [first, last]
    //! value constraints give pairwise-disjoint sets
    bool  valuesDisjoint(NodeIterator start, NodeIterator last);  /* TO DO */


    DictNode() {};                    /*< don't allow uninitialized node */
    Nodes            m_children;      /*< collection of child nodes */
    DictField       *m_field;         /*< id field for this node */
    DictNode        *m_parent;        /*< parent dictnode */
    DictConstraints *m_parConstraints;/*< Constraints on values of parent */
    DictConstraints *m_myConstraints; /*< Constraints on values for this node*/

  };       // end DictNode class definition
}          // end namespace xmlUtil
#endif

// $Header$
#ifndef XMLUTIL_DICTNODE_H
#define XMLUTIL_DICTNODE_H

#include "dom/DOM_Element.hpp"
#include "xmlUtil/id/DictObject.h"
#include "xmlUtil/id/DictField.h"
#include "xmlUtil/id/DictConstraints.h"
#include "xmlUtil/id/NamedId.h"

#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>

namespace xmlUtil {
  //! DictNode represents a node in an identifier dictionary.
  /*! It could be embedded in IdDict, but it's too big and complicated
      for that to be very attractive.  Instead, make virtually everything
      private and let IdDict be a friend.
  */
  class DictFieldMan;

  class DictNode : public DictObject {
  public:
    //! Construct a node from its XML representation (probably
    //! recursively construct children as well)
    DictNode(DOM_Element elt, DictNode* parent, DictFieldMan *fieldMan); 
    /* IN PROGRESS */
    ~DictNode();  // similarly, delete all children recursively /* TO DO */
    const DictField& getField() const {return *m_field;};

    // Following should be changed or add new function which
    //  returns local constraints if there are any,
    //  else if field has associated constraints, return pointer to that
    //  else return 0
    const DictConstraints& getValueConstraints() const { 
      return *m_myConstraints;};
    const DictConstraints& getParentConstraints() const {
      return *m_parConstraints;};
    bool accept(DictVisitor* vis);

    //! Check that collection of children define appropriately disjoint
    //! collections of Identifiers.  In the process may do some sorting
    //! of child nodes, so these functions are not const.  
    bool  consistentChildren();          

    //! Verify that parent constraints, if any, are consistent with
    //! other restrictions on values parent node may take on.
    bool  consistentParent();

    //! Verify that constraints (if any) on node's own value are
    //! consistent with constraints on field.
    bool  consistentValues();

    // Maybe also need function to check that field ref points
    // to actual field definition?

    // Perhaps need to add validity functions here to check
    //    parent constraints, if any, make sense
    //    own constraints, if any, make sense
    bool  allowed(const unsigned value) const;  
    //! Are the given child field name, child field value, and 
    //! value this node compatible with defn of node and children?
    bool  allowedChild(std::string childField, unsigned childValue, 
                       unsigned myValue) const;
    //! Are the given child field value, and value this node compatible 
    //! with defn of node and children?
    bool  allowedChild(unsigned childValue, unsigned myValue) const;

    //! Given an identifier, see if it is "allowed" by this node and
    //! a sequence of its descendents.  If there is a non-null
    //! NamedId argument, build up the associated NamedId
    bool allowIdentifier(const Identifier& id, NamedId* named = 0);

    //! return true iff there is a path starting with current node
    //! which could produce the NamedId \a nId
    bool allowNamedId(const NamedId& nId);

    //! return true iff there is a path starting with current node
    //! such that node field names match those in \a seq 
    bool allowNameSeq(const NameSeq& seq) const;

    //! Exception class
    class No_Assignment {};

  private:
    friend class IdDict;
    friend class IdConverter;   //maybe
    friend struct POrder;

    // Might want to use some other form of collection for 
    // child DictNodes
    typedef std::vector<DictNode*> Nodes;
    typedef std::vector<DictNode*>::iterator NodeIterator;
    typedef std::vector<DictNode*>::const_iterator ConstNodeIterator;

    // Nested class allowing sorting of child nodes by m_minValue of
    // parent constraint
    struct POrder {
      bool operator()(DictNode* a, DictNode* b) {
        return ((a->getParentConstraints()).getMin() <
                (b->getParentConstraints()).getMin()   ); }
    }; 

    const Nodes&     getChildren() const {return m_children;}

    /// The following may be called by IdDict during, e.g., validation
    static void   setOstream(std::ostream* out) {m_err = out;}

    // More constructor, for when dictionary is built dynamically
    // rather than from xml file
    /*    DictNode(DictField& field);                             TO DO
    DictNode(DictField& field, DictConstraints& parConstraints);  TO DO
    DictNode(DictField& field, DictConstraints& parConstraints,   TO DO
             DictConstraints& myConstraints);                     TO DO */

    //! Is the particular child node compatible with the given
    //! childValue and value for this?
    /*! Assume we've already checked that this node may take
        on myValue; still need to check that for the child indicated
        myValue is allowed.
    */
    bool  allowedChild(const DictNode* const thisChild, unsigned childValue,
                       unsigned myValue) const;


    //! Called by version without iterators
    bool allowIdentifier(Identifier::const_iterator idIt, 
                         Identifier::const_iterator end,
                         NamedId* named=0);

    //! Called by version with simple NamedId rather than iterators
    bool allowNamedId(NamedId::FieldIt nIdIt,
                      NamedId::FieldIt end);


    //! called by version with simple NameSeq arg rather than iterators
    bool allowNameSeq(NameSeq::const_iterator seqIt,
                      NameSeq::const_iterator end) const;

    bool  addChild(DictNode* child);


    //! Check that for nodes in the interval demarcated by [first, last]
    //! value constraints give pairwise-disjoint sets
    bool  valuesDisjoint(ConstNodeIterator start, 
                         ConstNodeIterator last); 


    //! Copy constructor copies all children, grandchildren, etc.
    DictNode(const DictNode& toCopy);

    // Assignment of nodes is not allowed
    DictNode& operator=(const DictNode&); 
    void deepCopy(const DictNode& toCopy);

    /** don't allow uninitialized node */
    DictNode() {};   

    Nodes            m_children;      /**< collection of child nodes */
    const DictField* m_field;         /**< id field for this node */
    const DictNode*  m_parent;        /**< parent dictnode */
    DictConstraints* m_parConstraints;/**< Constraints on values of parent */
    DictConstraints* m_myConstraints; /**< Constraints on values for this node*/
    static std::ostream* m_err;

  };       // end DictNode class definition
}          // end namespace xmlUtil
#endif

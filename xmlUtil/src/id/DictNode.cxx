// $Header$
#include "dom/DOM_Element.hpp"
#include "dom/DOMString.hpp"
#include "xml/Dom.h"
#include "xmlUtil/id/DictNode.h"
#include "xmlUtil/id/DictFieldMan.h"
#include <assert.h>

namespace xmlUtil {
  DictNode::DictNode(DOM_Element elt, DictNode *parent,
                     DictFieldMan *fieldMan) : 
    m_children(0), m_parent(parent), m_parConstraints(0), 
    m_myConstraints(0)   {
    // check element type is OK.  Shouldn't be a problem if
    // the xml file validates.  Note sense of assert is:
    //   abort program if assertion is *false*
    assert( (elt.getTagName().equals("dictNode")) ||
            (elt.getTagName().equals("dictRoot"))     ) ;

    // Field name stuff:
    std::string fName = 
      std::string(xml::Dom::transToChar(elt.getAttribute("fieldName")));
    
    m_field = fieldMan->find(fName.c_str());
    assert(m_field != 0);
    // Do something if it doesn't exist?? Shouldn't be terribly
    // necessary by ID/IDREF mechanism

    DOM_Element child = xml::Dom::getFirstChildElement(elt);
    // Could be a leaf with no children
    if (child == DOM_Element()) return;

    if ((child.getTagName()).equals("pValues")) { // parent constraints
      DOM_Element gChild = xml::Dom::getFirstChildElement(child);
      m_parConstraints = new DictConstraints(gChild);
      
      // check for compatibility
      DictConstraints* fromParent = 0;
      if (parent->m_myConstraints) fromParent = (parent->m_myConstraints);
      else {
        fromParent = parent->m_field->getConstraints();
      }
      if (fromParent) {
        bool allowed = (fromParent->allowed(m_parConstraints));
        if (!allowed) {
          assert(allowed);
        }
      }
      child = xml::Dom::getSiblingElement(child);
    }
    
    // Next up we might have our own constraint
    if (!(child.getTagName().equals("dictNode")) ) {
      m_myConstraints = new DictConstraints(child);
      
      // Check it's consistent with constraints on our own field,
      // if any
      DictConstraints* fCon = m_field->getConstraints();
      if (fCon != 0) {
        assert(fCon->allowed(m_myConstraints));
      }
      child = xml::Dom::getSiblingElement(child);
    }
    
    // All that's left are child nodes, if any
    
    while (child != DOM_Element()) {
      DictNode *nextChild = new DictNode(child, this, fieldMan);
      
      m_children.push_back(nextChild);
      child = xml::Dom::getSiblingElement(child);
    }
    // Finally, how about checking for consistency among children?
    if (!consistentChildren()) {
      // make a fuss
    }
  }
  
  // If 0 or 1 child, we're done.  Else 
  // *no* children have parent constraints
  //        or
  // *all* children must have parent constraints, 
  // and they must define sets
  // such that for any pair, the sets are disjoint or
  // they are identical.  Now collect all nodes with identical
  // parent constraints and check their "personal" value constraints.
  // They must be disjoint.
  bool DictNode::consistentChildren() {
    if (m_children.size() < 2) return true; // most common case
    
    // Either all children must have parent-value constraints or
    // none can.
    ConstNodeIterator it = m_children.begin();
    bool childMode = ((*it)->m_parConstraints != 0);
    ++it;
    for (; it < m_children.end(); ++it) {
      if (((*it)->m_parConstraints != 0) != childMode) return false;
    }
    if (!childMode) { // no parent constraints
      return valuesDisjoint(m_children.begin(), m_children.end());
    }
    
    // Otherwise sort nodes by minValue of parent constraint
    POrder parOrder;
    std::sort(m_children.begin(), m_children.end(), parOrder);
    
    // Check that parent constraints are disjoint or identical
    for (it = m_children.begin(); (it + 1) != m_children.end(); ++it) {
      DictConstraints *par1 = (*it)->m_parConstraints;
      DictConstraints *par2 = (*(it + 1))->m_parConstraints;
      
      if (  (!par1->equals(*par2)) &&
            (!par1->disjoint(*par2))  ) return false;
    }
    
    // Finally need to check that child values are disjoint
    // within subcollection of nodes with same parent constraints
    it = m_children.begin();
    while (it  != m_children.end() ) {
      ConstNodeIterator follow = it + 1;
      unsigned   ourMin = (*it)->m_parConstraints->getMin();
      
      // If we're ever starting a new collection to check at the
      // last child, we're done
      if (follow == m_children.end()) return true;
      while (follow != m_children.end()) {
        DictConstraints *par2 = (*follow)->m_parConstraints;
        
        if (par2->getMin() != ourMin) {
          --follow;
          break;
        }
        
      }
      // If we got here, it's time to check out a collection of children
      bool ok = valuesDisjoint(it, follow);
      if (!ok) return false;
      
      // Otherwise see if there is another bunch left
      it = ++follow;
    }
    return true;
  }

  bool DictNode::consistentParent() {
    if (!m_parConstraints) return true;

    if (DictConstraints *parConstraints = m_parent->m_myConstraints) {
      return parConstraints->allowed(*m_parConstraints);
    }
    else if (DictConstraints *parConstraints = 
             m_parent->m_field->getConstraints()) {
      return parConstraints->allowed(*m_parConstraints);
    }
    else return true;
  }

  bool DictNode::consistentValues() {
    if (!m_myConstraints) return true;

    DictConstraints *fieldConstraints = m_field->getConstraints();
    
    if (!fieldConstraints) return true;

    return fieldConstraints->allowed(*m_myConstraints);
  }
  
  bool DictNode::valuesDisjoint(ConstNodeIterator start, 
                                ConstNodeIterator last) {
    if (*start == *last) return true;
    
    ConstNodeIterator pCurrentNode = start;
    std::set<unsigned> values;
    (*start)->m_myConstraints->insertValues(values);

    while (pCurrentNode != last) {

      ++pCurrentNode;
      
      // Make a new set containing values for current node
      std::set<unsigned> curValues;
      (*pCurrentNode)->m_myConstraints->insertValues(curValues);

      // See if it intersects accumulated values from previous nodes
      std::vector<unsigned> intersect;
      std::set_intersection(values.begin(), values.end(), curValues.begin(),
                            curValues.end(), intersect.begin());
      if (intersect.size() > 0) return false;
      
      // If not, continue to accumulate
      (*pCurrentNode)->m_myConstraints->insertValues(values);

    } 
    return true;
    
  }
  
  bool DictNode::allowed(const unsigned value) const {
    if (!m_myConstraints) {
      return m_field->allowed(value);
    }
    else {
      return m_myConstraints->allowed(value);
    }
  }

  bool  DictNode::allowedChild(std::string childField, unsigned childValue, 
                               unsigned myValue) const {
    if (!allowed(myValue)) return false;

    for (ConstNodeIterator it = m_children.begin(); it != m_children.end();
         ++it) {
      if (  ((*it)->getField()).getName() == childField) {
        return allowedChild(*it, childValue, myValue);
      }
    }
    // Never found child with correct field name
    return false;
  }
          
  bool  DictNode::allowedChild(unsigned childValue, 
                               unsigned myValue) const {

    if (!allowed(myValue)) return false;
    for (ConstNodeIterator it = m_children.begin(); it != m_children.end();
         ++it) {
      if  (allowedChild(*it, childValue, myValue)) return true;
    }
    return false;
  }

  bool DictNode::allowedChild(const DictNode* const child, unsigned childValue,
                              unsigned myValue) const {
    if (child->m_parConstraints) {
      if (!(child->m_parConstraints)->allowed(myValue)) return false;
    }
    return child->allowed(childValue);
  }
  bool DictNode::allowIdentifier(Identifier::iterator idIt, 
                                 Identifier::const_iterator end) {
    if (!allowed(*idIt)) return false;
    Identifier::iterator tmp = idIt;
    ++tmp;
    if (tmp == end) return true;

    //look for children allowing our value for parent node
    for (ConstNodeIterator it = m_children.begin(); it != m_children.end();
         it++) {
      DictNode *child = *it;
      if ( (child->m_parConstraints == 0) ||
           (child->m_parConstraints->allowed(*idIt)) ) {
        if (child->allowIdentifier(++idIt, end)) return true;
      }
    }
    return false;
  }
    

  bool DictNode::allowIdentifier(const Identifier& id) {
    Identifier::iterator idIt = id.begin();
    return allowIdentifier(idIt, id.end());
  }

  bool DictNode::allowNamedId(const NamedId& nId) {
    NamedId::FieldIt nIdIt = nId.m_fields->begin();
    return allowNamedId(nIdIt, nId.m_fields->end());
  }

  bool DictNode::allowNamedId(NamedId::FieldIt nIdIt,
                              NamedId::FieldIt end) {
    // Check information for first field matches

    // These aren't supposed to happen
    if (nIdIt == end) return true;
    if (!(*nIdIt)) return false; 

    if ((*nIdIt)->name != m_field->getName()) return false;
    if (!allowed((*nIdIt)->value) ) return false;

    NamedId::FieldIt tmp = nIdIt;
    ++tmp;
    if (tmp == end) return true; //done

    // Now search children for match of next field
    for (ConstNodeIterator it = m_children.begin(); it != m_children.end();
         it++) {
      DictNode *child = *it;
      if ( (child->m_parConstraints == 0) ||
           (child->m_parConstraints->allowed((*nIdIt)->value) ) ) {
        if (child->allowNamedId(++nIdIt, end) ) return true;
      }
    }
    return false;
  }

  bool  DictNode::addChild(DictNode* child) {
    if ((m_myConstraints != 0) & (child->m_parConstraints != 0)) {
      if (!(m_myConstraints->allowed(child->m_parConstraints))) return false;
    }
    m_children.push_back(child);
    return true;
  }


  // Support visitors.
  // Should there be a non-recursive version?
  // Should there be an option to visit associated field, constraints?
  bool DictNode::accept(DictVisitor* vis) {
    // Call back visitor for child nodes, then self

    for (ConstNodeIterator it = m_children.begin(); it != m_children.end();
         it++) {
      bool keepGoing = (*it)->accept(vis);
      if (!keepGoing) return false;
    } 

    return  vis->visitNode(this);
  }

  DictNode::DictNode(const DictNode& toCopy) : DictObject(toCopy) {
    deepCopy(toCopy);
  }

  DictNode& DictNode::operator=(const DictNode& d) {
    throw No_Assignment();
  }

  void DictNode::deepCopy(const DictNode& d) {
    // Don't copy parent.  Nodes created with copy constructor
    // are parentless.  Use addChild after creation.
    // Also have to get parent constraints some other way
    m_parent = 0;
    m_parConstraints = 0;
    if (d.m_field) {
      m_field = new DictField(*(d.m_field));
    }
    else m_field = 0;

    if (d.m_myConstraints) {
      m_myConstraints = new DictConstraints(*(d.m_myConstraints));
    }
    else m_myConstraints = 0;

    // Now for child nodes
    for (ConstNodeIterator it = d.m_children.begin(); 
         it != d.m_children.end(); ++it) {
      DictNode *newChild = new DictNode(**it);
      DictConstraints *oldParConstraints = (*it)->m_parConstraints;

      newChild->m_parConstraints = (oldParConstraints) ? 
        new DictConstraints(*oldParConstraints) : 0;
      // Need some exception handling here.  addChild returns 
      // true or false, depending on whether child actually was
      // added or not.
      addChild(newChild);
    }
  }

  DictNode::~DictNode() {
    if (m_parConstraints != 0) {
      delete m_parConstraints;
      m_parConstraints = 0;
    }
    if (m_myConstraints != 0) {
      delete m_myConstraints;
      m_myConstraints = 0;
    }
    for (NodeIterator it = m_children.begin(); 
         it != m_children.end(); ++it) {
      delete *it;
    }
  }
}                  // end namespace



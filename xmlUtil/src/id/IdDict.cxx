// $Header$

#include "dom/DOMString.hpp"
#include "dom/DOM_NodeList.hpp"
#include "xml/Dom.h"
#include "xmlUtil/id/IdDict.h"
#include "xmlUtil/id/DictFieldMan.h"
#include "xmlUtil/id/DictField.h"
#include "xmlUtil/id/DictNode.h"
#include "xmlUtil/id/DictValidVisitor.h"
#include <assert.h>
namespace xmlUtil {
  IdDict::IdDict(DOM_Element elt)  {
    // Check that element has the right tag name: idDict
    // Caller probably will have done this already
    assert(elt.getTagName().equals("idDict"));

    // Store information from attributes
    m_name = std::string(xml::Dom::getAttribute(elt, "name"));

    std::string intVal = xml::Dom::getAttribute(elt, "major");
    m_major = atoi(intVal.c_str());

    intVal = xml::Dom::getAttribute(elt, "minor");
    m_minor = atoi(intVal.c_str());

    intVal = xml::Dom::getAttribute(elt, "patch");
    m_patch = atoi(intVal.c_str());

    // Check number of children.  This is an upper bound
    // on number of fields we have to store
    int size = (elt.getChildNodes()).getLength();
    m_fieldMan = new DictFieldMan(size);
                
    DOM_Element fieldElt = xml::Dom::getFirstChildElement(elt);

    // Register all the fields
    while ((fieldElt.getNodeName()).equals("field")) {
      //      DOM_Element fieldElt = static_cast<DOM_Element&>(fieldNode);
      DictField* field = new DictField(fieldElt);
      m_fieldMan->signup(field);
      fieldElt = xml::Dom::getSiblingElement(fieldElt);
    }
    // Finally make the hierarchy of constraints on Identifiers
    m_root = new DictNode(fieldElt, 0, m_fieldMan);
  }

  IdDict::~IdDict() {
    /* First get rid of hierarchy of nodes, then get rid of
       our DictFieldMan object. The latter will delete all the
       fields
    */
    delete m_root;
    delete m_fieldMan;
  }

  bool IdDict::accept(DictVisitor *vis) {
    return accept(vis, 0xffffffff);
  }


  bool IdDict::accept(DictVisitor *vis, unsigned mask) {
    bool status = true;

    if (!vis->visitDict(this)) return false;

    if (mask & nodeHierarchy) 
      status = m_root->accept(vis);

    if (status && (mask & fieldManager)) {
      status &= m_fieldMan->accept(vis);
    }
    return status;
  }


  bool IdDict::isValid()  {
    DictValidVisitor visitor;

    accept(&visitor);
    return visitor.wasValid();
  }

  bool IdDict::idOk(const Identifier& id) const {
    return m_root->allowIdentifier(id);
  }

  bool IdDict::idOk(const NamedId& id) const {
    return m_root->allowNamedId(id);
  }

  bool IdDict::nameSeqOk(const NameSeq& seq) const {
    return m_root->allowNameSeq(seq);
  }

  // all the work is in DictNode::allowIdentifier(const Identifier&, NamedId*)
  NamedId *  IdDict::getNamedId(const Identifier& id) const {
    NamedId *named = new NamedId;
    m_root->allowIdentifier(id, named);
    return named;
  }

  bool IdDict::addChild(DictNode* parent, DictNode* newNode) {
    /* Started to write something to add field if not already
       registered, but whoever has created newNode should
       already have done this */
    return parent->addChild(newNode);
  }


}

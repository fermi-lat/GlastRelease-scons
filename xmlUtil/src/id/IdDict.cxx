// $Header$

//#include <xercesc/dom/DOMString.hpp>
//#include <xercesc/dom/DOM_NodeList.hpp>
#include "xmlBase/Dom.h"
#include "xmlUtil/id/IdDict.h"
#include "xmlUtil/id/DictFieldMan.h"
#include "xmlUtil/id/DictField.h"
#include "xmlUtil/id/DictNode.h"
#include "xmlUtil/id/DictValidVisitor.h"
#include <assert.h>
namespace xmlUtil {
  XERCES_CPP_NAMESPACE_USE
  IdDict::IdDict(DOMElement* elt)  {
    // Check that element has the right tag name: idDict
    // Caller probably will have done this already
    if (!xmlBase::Dom::checkTagName(elt, "idDict")) 
      throw xmlBase::DomException("Invalidly constructed id dictionary");

    //    assert(elt.getTagName().equals("idDict"));

    // Store information from attributes
    m_name = std::string(xmlBase::Dom::getAttribute(elt, "name"));

    try {
      m_major = xmlBase::Dom::getIntAttribute(elt, "major");
      m_minor = xmlBase::Dom::getIntAttribute(elt, "minor");
      m_patch = xmlBase::Dom::getIntAttribute(elt, "patch");
    }
    catch (xmlBase::WrongAttributeType ex) {
      std::cerr << "From xmlUtil::IdDict::IdDict  bad version specifier " 
                << std::endl << ex.getMsg() << std::endl;
      throw ex;
    }

    // Check number of field children.
    std::vector<DOMElement*> children;
    xmlBase::Dom::getChildrenByTagName(elt, "field", children);
    unsigned int nField = children.size();
    m_fieldMan = new DictFieldMan(nField);
                
    for (unsigned iField = 0; iField < nField; iField++) {
      // Register all the fields
      DictField* field = new DictField(children[iField]);
      m_fieldMan->signup(field);
    }
    //DictRoot immediately follows last field element
    DOMElement* dictRoot = xmlBase::Dom::getSiblingElement(children.back());
    // Finally make the hierarchy of constraints on Identifiers
    m_root = new DictNode(dictRoot, 0, m_fieldMan);
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

  NameSeq*   IdDict::getNameSeq(const Identifier& id) const {
    NameSeq *seq = new NameSeq;
    m_root->allowIdentifier(id, seq);
    return seq;
  }

  std::string IdDict::getNameSeqString(const Identifier& id) const {
    if (id.size() == 0) {
      std::string s = "";
      return s;
    }
    NameSeq* seq = getNameSeq(id);
    std::string s = nameSeqString(*seq);
    delete seq;
    return s;
  }

  bool IdDict::addChild(DictNode* parent, DictNode* newNode) {
    /* Started to write something to add field if not already
       registered, but whoever has created newNode should
       already have done this */
    return parent->addChild(newNode);
  }
}

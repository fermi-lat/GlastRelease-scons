// $Header$

#include <vector>
#include <string>
#include "xmlUtil/Substitute.h"
#include "xmlBase/Dom.h"
#include "xmlUtil/Arith.h"
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMDocument.hpp>
// #include <xercesc/dom/DOMTreeWalker.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>
#include <xercesc/dom/DOMAttr.hpp>
#include <xercesc/util/XMLString.hpp>
#include <vector>

namespace xmlUtil {
  XERCES_CPP_NAMESPACE_USE

  Substitute::Substitute(DOMDocument* doc, std::string suffix) :
    m_doc(doc), m_count(0), m_suffix(suffix) {
    m_suffixLen = suffix.size();
  }

  int Substitute::execute(DOMElement* top) {
    // Get rid of some inefficiency by means of the following imperfect
    // mechanism:
    //   An element may have an attribute called "substituted".  If it
    //   does, the attribute always has a default value of "false".
    //   The execute method looks for the attribute.  If it is there and
    //   has value "true", substitution has already been done, so just
    //   return.  If the attribute has value "false" or doesn't exist,
    //   go ahead and do the substitution.  Afterwards, if the attribute
    //   exists set its value to "true", but don't add any new attributes.
    //   Also as a matter of course always check the root element for
    //   a "substituted" attribute.
    //
    //   If an element and its descendent *both* have the substituted
    //   attribute, only the one passed to execute will have its substitute
    //   attribute examined and possibly changed.  This means it's possible
    //   for an element to have a "substituted" attribute value of "false"
    //   when it really should be "true", but never the other way around.
    //   
    DOMElement* docElt = m_doc->getDocumentElement();

    if (std::string("true") == 
        xmlBase::Dom::getAttribute(docElt,"substituted")) return 0;
    
    if (std::string("true") == xmlBase::Dom::getAttribute(top, "substituted"))
      return 0;

    // This is not explained anywhere in the xerces doc. that I
    // can see, but from the code in TreeWalkerImpl::acceptNode
    // this would seem to be the right value for whatToShow
    //    unsigned long whatToShow = 1 << (DOMNode::ELEMENT_NODE -1);

    std::vector<DOMElement*> allElts;
    xmlBase::Dom::getDescendantsByTagName(top, "*", allElts);

    unsigned int nElt = allElts.size();

    for (unsigned int iElt = 0; iElt < nElt; iElt++) {
      sub(allElts[iElt]);
    }
    XMLCh* xmlchSubs = XMLString::transcode("substituted");
    DOMAttr* subNode = top->getAttributeNode(xmlchSubs);
    XMLString::release(&xmlchSubs);
    if (subNode != 0) { 
      XMLCh* xmlchTrue = XMLString::transcode("true");
      subNode->setValue(xmlchTrue);
      XMLString::release(&xmlchTrue);
    }

    return m_count;
  }

  void Substitute::sub(DOMElement* elt) {
    using xmlBase::Dom;

    DOMNamedNodeMap* attMap = elt->getAttributes();
    int   nAtt = attMap->getLength();
    int   iAtt;
    std::vector<std::string> toProcess;
    std::vector<unsigned> toProcessPos;

    for (iAtt = 0; iAtt <nAtt; iAtt++) {
      DOMNode* att = attMap->item(iAtt);
      std::string attName = Dom::getNodeName(att);
      
      // Look for suffix
      unsigned pos = attName.find(m_suffix, attName.size() - m_suffixLen);
      if  (pos < attName.size())  { // found
        toProcess.push_back(attName);
        toProcessPos.push_back(pos);
      }
    }

    // Now process by name
    while (toProcess.size() > 0) {
      std::string oldAttName = toProcess.back();
      unsigned  pos = toProcessPos.back();
      toProcess.pop_back();
      toProcessPos.pop_back();
      std::string eltId = Dom::getAttribute(elt, oldAttName);
      DOMElement* constElt =  Dom::getElementById(m_doc, eltId);

      if (constElt == 0 ) { // shouldn't happen
        m_notFound++;
        continue;
      }
      // Build the new attribute
      std::string val = Dom::getAttribute(constElt, "value");

      // If the element referred to is a <prim>, we should be done
      // because it always has a value.  If it's a <const>
      // we may have to evaluate it.
      if (!val.size()) { // try evaluating it
        Arith curArith(constElt);
        /*double evalValue = */
        curArith.evaluate();
        curArith.saveValue();
        val = Dom::getAttribute(constElt, "value");
      }

      std::string newAttName = std::string(oldAttName, 0, pos);

      Dom::addAttribute(elt, newAttName, val);
      XMLCh* xmlchOldAttName = XMLString::transcode(oldAttName.c_str());
      elt->removeAttribute(xmlchOldAttName);
      XMLString::release(&xmlchOldAttName);
      m_count++;
    }

  }
}

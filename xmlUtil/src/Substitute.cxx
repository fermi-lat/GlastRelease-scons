// $Header$

#include <vector>
#include <string>
#include "xmlUtil/Substitute.h"
#include "xml/Dom.h"
#include "xmlUtil/Arith.h"
#include <xercesc/dom/DOM_Element.hpp>
#include <xercesc/dom/DOM_TreeWalker.hpp>
#include <xercesc/dom/DOM_NamedNodeMap.hpp>
#include <vector>

namespace xmlUtil {

  Substitute::Substitute(DOM_Document doc, std::string suffix) :
    m_doc(doc), m_count(0), m_suffix(suffix) {
    m_suffixLen = suffix.size();
  }

  int Substitute::execute(DOM_Element top) {
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
    DOM_Element docElt = m_doc.getDocumentElement();

    if (std::string("true") == xml::Dom::getAttribute(docElt, "substituted"))
      return 0;
    
    if (std::string("true") == xml::Dom::getAttribute(top, "substituted"))
      return 0;

    DOM_Node curNode = top;

    // This is not explained anywhere in the xerces doc. that I
    // can see, but from the code in TreeWalkerImpl::acceptNode
    // this would seem to be the right value for whatToShow
    unsigned long whatToShow = 1 << (DOM_Node::ELEMENT_NODE -1);

    DOM_TreeWalker walker = m_doc.createTreeWalker(top, whatToShow, 0, 0);
    m_count = 0;
    while (curNode != DOM_Node() ) {
      sub(static_cast<DOM_Element &> (curNode));
      curNode = walker.nextNode();
    }
    DOM_Attr subNode = top.getAttributeNode("substituted");
    if (subNode != DOM_Attr()) { 
      subNode.setValue("true");
    }

    return m_count;
  }

  void Substitute::sub(DOM_Element elt) {
    using xml::Dom;

    DOM_NamedNodeMap attMap = elt.getAttributes();
    int   nAtt = attMap.getLength();
    int   iAtt;
    std::vector<std::string> toProcess;
    std::vector<unsigned> toProcessPos;

    for (iAtt = 0; iAtt <nAtt; iAtt++) {
      DOM_Node att = attMap.item(iAtt);
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
      DomElement constElt = 
        Dom::getElementById(m_doc, Dom::getAttribute(elt, oldAttName));

      if (constElt == DomElement() ) { // shouldn't happen
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
        double evalValue = curArith.evaluate();
        curArith.saveValue();
        val = Dom::getAttribute(constElt, "value");
      }

      std::string newAttName = std::string(oldAttName, 0, pos);

      Dom::addAttribute(elt, newAttName, val);
      elt.removeAttribute(oldAttName.c_str());
      m_count++;
    }

  }
}

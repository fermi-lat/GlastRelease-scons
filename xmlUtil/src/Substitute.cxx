// $Header$

#include <vector>
#include <string>
#include "xmlUtil/Substitute.h"
#include "xml/Dom.h"
#include "xmlUtil/Arith.h"
#include "dom/DOM_Element.hpp"
#include "dom/DOM_TreeWalker.hpp"
#include "dom/DOM_NamedNodeMap.hpp"
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
    if ((docElt.getAttribute("substituted")).equals(DOMString("true"))) 
      return 0;
    
    if ((top.getAttribute("substituted")).equals(DOMString("true"))) return 0;
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
    DOM_NamedNodeMap attMap = elt.getAttributes();
    int   nAtt = attMap.getLength();
    int   iAtt;
    std::vector<DOMString> toProcess;
    std::vector<unsigned> toProcessPos;
    //    std::vector<DOM_Node> toDelete;

    for (iAtt = 0; iAtt <nAtt; iAtt++) {
      DOM_Node att = attMap.item(iAtt);
      //      if (att != DOM_Node()) { // probably not necessary
      DOMString  attString = att.getNodeName();
      std::string attName = std::string(xml::Dom::transToChar(attString));
      
        // Look for suffix
      unsigned pos = attName.find(m_suffix, attName.size() - m_suffixLen);
      if  (pos < attName.size())  { // found
        toProcess.push_back(attString);
        toProcessPos.push_back(pos);
      }
    }

    // Now process by name
    while (toProcess.size() > 0) {
      DOMString oldAttString = toProcess.back();
      std::string oldAttName = 
        std::string(xml::Dom::transToChar(oldAttString));
      unsigned  pos = toProcessPos.back();
      toProcess.pop_back();
      toProcessPos.pop_back();
      DOM_Element constElt =
        m_doc.getElementById(elt.getAttribute(oldAttString));
      if (constElt == DOM_Element() ) { // shouldn't happen
        m_notFound++;
        continue;
      }
      // Build the new attribute
      DOMString val = constElt.getAttribute(DOMString("value"));

      // If the element referred to is a <prim>, we should be done
      // because it always has a value.  If it's a <const>
      // we may have to evaluate it.
      if (val == DOMString()) { // try evaluating it
        Arith curArith(constElt);
        double evalValue = curArith.evaluate();
        curArith.saveValue();
        val = constElt.getAttribute(DOMString("value"));
      }

      std::string newAttName = oldAttName.erase(pos, oldAttName.size()  );
      elt.setAttribute(DOMString(newAttName.c_str()), val);
      elt.removeAttribute(oldAttString);
      m_count++;
    }

    // Now delete each saved attribute node 
    // while (toDelete.size() > 0) {
    //      DOM_Node att = toDelete.back();
    //      toDelete.pop_back();
    //      elt.removeAttributeNode(static_cast <DOM_Attr &>(att));
    //    }
  }
}

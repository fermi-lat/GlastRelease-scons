// $Header$

#include <string>
#include "xmlUtil/Arith.h"
#include "xml/Dom.h"
#include <dom/DOM_Document.hpp>

namespace xmlUtil {

  inline DOM_Element firstEltChild(DOM_Element& elt) {
    DOM_Node curNode = elt.getFirstChild();
    if (curNode == DOM_Node()) return DOM_Element();
    while (curNode.getNodeType() != DOM_Node::ELEMENT_NODE) {
      curNode = curNode.getNextSibling();
      if (curNode == DOM_Node()) return DOM_Element();
    }
    return DOM_Element(static_cast<DOM_Element &>(curNode));
  }
  
  inline  DOM_Element nextEltSibling(DOM_Element& elt) {
    DOM_Node curNode = elt.getNextSibling();
    if (curNode == DOM_Node()) return DOM_Element();
    while (curNode.getNodeType() != DOM_Node::ELEMENT_NODE) {
      curNode = curNode.getNextSibling();
      if (curNode == DOM_Node()) return DOM_Element();
    }
    return DOM_Element(static_cast<DOM_Element &>(curNode));
  }

  // Define statics here.
  // Rather than continually creating and destroying DOMStrings
  // needed for comparison, initialize them once and for all in 
  // init()
  Arith::ptrDOMString * Arith::typeNames;
  Arith::ptrDOMString Arith::valString;   // a DOMString containing "value"
  Arith::ptrDOMString Arith::refToString;
  Arith::ptrDOMString Arith::constString;
  Arith::ptrDOMString Arith::referString;
  Arith::ptrDOMString Arith::addString;
  Arith::ptrDOMString Arith::minusString;
  Arith::ptrDOMString Arith::uminusString;
  Arith::ptrDOMString Arith::mulString;
  Arith::ptrDOMString Arith::quoString;
  Arith::ptrDOMString Arith::maxString;
  Arith::ptrDOMString Arith::notesString;
  Arith::ptrDOMString Arith::lengthString;
  Arith::ptrDOMString Arith::cmString;
  Arith::ptrDOMString Arith::mString;
  
  Arith::Arith(const DOM_Element elt) {
    int i = 0;
    bool notFound = true;
    DOMString tagName = elt.getTagName();
    m_elt = elt;
    m_evaluated = false;
    m_tag = -1;
    if (typeNames == 0) init();

    while ((notFound) && (i < ETAG_n)) {
      DOMString *p = typeNames[i];
      DOMString tmp = DOMString(*p);
      if (tagName.equals(tmp) ) {
        m_tag = i;
        notFound = false;
      }
      i++;
    }
  }

  double Arith::evaluate() {
    if (!m_evaluated) { // evaluate
      DOM_Element curElt;
      switch(m_tag) {
      case ETAG_const: {
        DOMString val = m_elt.getAttribute(*valString);
        if (!(val.equals(DOMString()) ) ) { //got something 
          // Had better convert properly to a double
          m_number = atof(xml::Dom::transToChar(val));
          m_number *= getUnitsScale(m_elt);
        }
        else  { // must have a single operator child or refer child
          // either of which can be evaluated, optionally preceded
          // by a <notes> child, which we ignore
          curElt = firstEltChild(m_elt);
          if (curElt.getTagName().equals(*notesString)) { // move on
            DOM_Node next = curElt.getNextSibling();
            curElt = static_cast <DOM_Element &>(next);
          }
          m_number = Arith(curElt).evaluate();
        }
        break;
      }
      case ETAG_refer: {
        // Find the element pointed to and evaluate it
        DOMString ref = m_elt.getAttribute(*refToString);
        curElt = (m_elt.getOwnerDocument()).getElementById(ref);
        if (curElt.getTagName().equals("prim")) {
          m_number = 
            atof(xml::Dom::transToChar(curElt.getAttribute(*valString)));
          m_number *= getUnitsScale(curElt);
        }
        else {
          m_number = Arith(curElt).evaluate();
        }
        break;
      }
      case ETAG_uminus: {
        // Have a single child
        curElt = firstEltChild(m_elt);
        m_number = -(Arith(curElt).evaluate());
        break;
      }
      // Do the work elsewhere for the remaining operators
      case ETAG_add: {
        m_number = add();
        break;
      }
      
      case ETAG_mul: {
        m_number = mul();
        break;
      }
      
      case ETAG_minus: {
        m_number = minus();
        break;
      }
      case ETAG_quo: {
        m_number = quo();
        break;
      }
      
      case ETAG_max: {
        m_number = myMax();
        break;
      }
      default: {
        return 0;
      }
      }
      m_evaluated = true;
    }
    return m_number;
  }

  double Arith::add() {
    double  ans = 0.0;
    DOM_Element elt = firstEltChild(m_elt);
    ans += Arith(elt).evaluate();
    elt = nextEltSibling(elt);
    while (elt != DOM_Element()) {
      ans += Arith(elt).evaluate();
      elt = nextEltSibling(elt);
    }
    return ans;
  }

  double Arith::mul() {
    double  ans = 1.0;
    DOM_Element elt = firstEltChild(m_elt);
    ans *= Arith(elt).evaluate();
    elt = nextEltSibling(elt);
    while (elt != DOM_Element()) {
      ans *= Arith(elt).evaluate();
      elt = nextEltSibling(elt);
    }
    return ans;
  }

  double Arith::myMax() {
    DOM_Element elt = firstEltChild(m_elt);
    double  ans = Arith(elt).evaluate();
    elt = nextEltSibling(elt);
    while (elt != DOM_Element()) {
      double newAns = Arith(elt).evaluate();
      if (newAns > ans) ans = newAns;
      elt = nextEltSibling(elt);
    }
    return ans;
  }
  double Arith::minus() {
    double  ans;
    DOM_Element elt = firstEltChild(m_elt);
    ans = Arith(elt).evaluate();
    elt = nextEltSibling(elt);
    ans -= Arith(elt).evaluate();
    return ans;
  }

  double Arith::quo() {
    double  ans;
    double  divisor;
    DOM_Element elt = firstEltChild(m_elt);
    ans = Arith(elt).evaluate();
    elt = nextEltSibling(elt);
    divisor = Arith(elt).evaluate();  // check for 0?
    ans /= divisor;
    return ans;
  }

  // Store value back in DOM (evaluating first if necessary)
  void Arith::saveValue() {
    // Only makes sense for a const
    if (m_tag != ETAG_const) return;

    if (!m_evaluated) evaluate();

    // Check if we're supposed to be an int.  If so, coerce
    // m_number to be nearby int in case of round-off error 
    if (DOMString("int").equals(m_elt.getAttribute("type"))) {
      // If we're not already a perfect int, attempt to fix
      // so that we round the right way.
      // Otherwise, leave well enough alone
      long int intValue = m_number;
      double   intified = intValue;
      if (intified != m_number) {
        double fixup = 0.5;
        if (m_number < 0.0) fixup = -fixup;
        long int  intValue = (m_number + fixup);  
        m_number = intValue;
      }
    }
    // Answer will always be returned in mm for length-type
    // constants, so, just in case, this attribute should be
    // set to "mm".  If we're not a length-type constant,
    // setting this attribute is harmless.

    xml::Dom::addAttribute(m_elt, *valString, m_number);
    xml::Dom::addAttribute(m_elt, std::string("unitLength"), 
                           std::string("mm"));
  }

  double Arith::getUnitsScale(const DOM_Element& elt) {
    if (!(elt.getAttribute("uType")).equals(*lengthString) ) {
      return 1.0;
    }

    DOMString unitLength = elt.getAttribute("unitLength");
    if (unitLength.equals(*cmString)) return 10.0;
    else if (unitLength.equals(*mString)) return 1000.0;
    else return 1.0;
  }

  void Arith::init() {
    constString = new DOMString("const");
    referString = new DOMString("refer");
    addString = new DOMString("add");
    minusString = new DOMString("minus");
    mulString = new DOMString("mul");
    quoString = new DOMString("quo");
    uminusString = new DOMString("uminus");
    maxString = new DOMString("max");
    notesString = new DOMString("notes");
    lengthString = new DOMString("length");
    cmString = new DOMString("cm");
    mString = new DOMString("m");

    valString = new DOMString("value");
    refToString = new DOMString("refTo");


    typeNames = new ptrDOMString[ETAG_n];
    typeNames[ETAG_const] = constString;
    typeNames[ETAG_refer] = referString;
    typeNames[ETAG_add] = addString;
    typeNames[ETAG_minus] = minusString;
    typeNames[ETAG_mul] = mulString;
    typeNames[ETAG_quo] = quoString;
    typeNames[ETAG_uminus] = uminusString;
    typeNames[ETAG_max] = maxString;
  }
}

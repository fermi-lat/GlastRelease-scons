// $Header$

// #include <string>
#include "xmlUtil/Arith.h"
#include "xml/Dom.h"
#include <xercesc/dom/DOM_Document.hpp>
#include <string.h>  // for comparison of old-style c strings.

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
  Arith::ptrString * Arith::typeNames;
  Arith::ptrString Arith::valString = "value"; 
  Arith::ptrString Arith::refToString = "refTo";
  Arith::ptrString Arith::constString = "const";
  Arith::ptrString Arith::referString = "refer";
  Arith::ptrString Arith::addString = "add";
  Arith::ptrString Arith::minusString = "minus";
  Arith::ptrString Arith::uminusString = "uminus";
  Arith::ptrString Arith::mulString = "mul";
  Arith::ptrString Arith::quoString = "quo";
  Arith::ptrString Arith::maxString = "max";
  Arith::ptrString Arith::notesString = "notes";
  Arith::ptrString Arith::lengthString = "length" ;
  Arith::ptrString Arith::cmString = "cm";
  Arith::ptrString Arith::mString = "m";
  Arith::ptrString Arith::halfString = "half";
  
  Arith::Arith(const DOM_Element elt) {
    int i = 0;
    bool notFound = true;
    //    std::string tagName = xml::Dom::getTagName(elt);
    const char* tagName = (xml::Dom::getTagName(elt)).c_str();
    m_elt = elt;
    m_evaluated = false;
    m_tag = -1;
    if (typeNames == 0) init();

    while ((notFound) && (i < ETAG_n)) {
      ptrString p = typeNames[i];
      /*      std::string tmp = std::string(*p); */
      /*      if (tagName.equals(tmp) ) { */
      if (!strcmp(p, tagName)) {
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
        if (xml::Dom::hasAttribute(m_elt, valString)) {
          try {
            m_number = 
              xml::Dom::getDoubleAttribute(m_elt, std::string(valString) );
          }
          catch(xml::WrongAttributeType ex) {
            std::cerr << std::endl << ex.getMsg() << std::endl;
            throw ex;
          }
          m_number *= getUnitsScale(m_elt);
        }
        else  { // must have a single operator child or refer child
          // either of which can be evaluated, optionally preceded
          // by a <notes> child, which we ignore
          curElt = firstEltChild(m_elt);

          if (std::string(notesString) == xml::Dom::getTagName(curElt)) {
            DOM_Node next = curElt.getNextSibling();
            curElt = static_cast <DOM_Element &>(next);
          }
          m_number = Arith(curElt).evaluate();
        }
        break;
      }
      case ETAG_refer: {
        // Find the element pointed to and evaluate it
        std::string ref = xml::Dom::getAttribute(m_elt, refToString);
        curElt = (m_elt.getOwnerDocument()).getElementById(ref.c_str());

        if (std::string("prim") == xml::Dom::getTagName(curElt) ) {

          try {
            m_number = 
              xml::Dom::getDoubleAttribute(curElt, std::string(valString) );
          }
          catch(xml::WrongAttributeType ex) {
            std::cerr << std::endl << ex.getMsg() << std::endl;
            throw ex;
          }

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
      case ETAG_half: {
        // Have a single child
        curElt = firstEltChild(m_elt);
        m_number = 0.5 * (Arith(curElt).evaluate());
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

    // Used to check if we were supposed to be an int and if so, coerce
    // m_number to be nearby int in case of round-off error, but this
    // shouldn't actually happen -- in no case do we do arithmetic
    // and attempt to produce an int result.  Code has been excised.

    if (xml::Dom::getAttribute(m_elt, "type") == std::string("int") ) {
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

    xml::Dom::addAttribute(m_elt, std::string(valString), m_number);
    xml::Dom::addAttribute(m_elt, std::string("unitLength"), 
                           std::string("mm"));
  }

  double Arith::getUnitsScale(const DOM_Element& elt) {
    if (xml::Dom::getAttribute(elt, "uType") != std::string(lengthString)) {
      return 1.0;
    }

    
    std::string unitLength = xml::Dom::getAttribute(elt, "unitLength");
    if (unitLength == std::string(cmString)) return 10.0;
    else if (unitLength == std::string(mString)) return 1000.0;
    else return 1.0;
  }

  void Arith::init() {
    typeNames = new ptrString[ETAG_n];
    typeNames[ETAG_const] = constString;
    typeNames[ETAG_refer] = referString;
    typeNames[ETAG_add] = addString;
    typeNames[ETAG_minus] = minusString;
    typeNames[ETAG_mul] = mulString;
    typeNames[ETAG_quo] = quoString;
    typeNames[ETAG_uminus] = uminusString;
    typeNames[ETAG_max] = maxString;
    typeNames[ETAG_half] = halfString;
  }
}

// $Header$

// #include <string>
#include "xmlBase/Dom.h"
#include "xmlUtil/Arith.h"
#include <xercesc/dom/DOMDocument.hpp>
// #include <xercesc/dom/DOMElement.hpp>

    
    /* keep some static strings so we don't have to build them over and over */

namespace xmlUtil {
  //  using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
  //  using XERCES_CPP_NAMESPACE_QUALIFIER DOMNode;
  XERCES_CPP_NAMESPACE_USE

  Arith::ptrString * Arith::typeNames = 0;
  Arith::ptrString Arith::valString = 0; // pointer to std::string of  "value"
  Arith::ptrString Arith::refToString = 0; // ..and so forth
  Arith::ptrString Arith::notesString = 0;
  Arith::ptrString Arith::lengthString = 0;
  Arith::ptrString Arith::cmString = 0;
  Arith::ptrString Arith::mString = 0;

 void Arith::init() {
    typeNames = new ptrString[ETAG_n];

    typeNames[ETAG_const] = new std::string("const");
    typeNames[ETAG_refer] = new std::string("refer");
    typeNames[ETAG_add] = new std::string("add");
    typeNames[ETAG_minus] = new std::string("minus");
    typeNames[ETAG_mul] = new std::string("mul");
    typeNames[ETAG_quo] = new std::string("quo");
    typeNames[ETAG_uminus] = new std::string("uminus");
    typeNames[ETAG_max] = new std::string("max");
    typeNames[ETAG_half] = new std::string("half");

    notesString = new (std::string) ("notes");
    lengthString = new std::string("length");
    refToString = new std::string("refTo");
    cmString = new std::string("cm");
    mString = new std::string("m");
    valString = new std::string("value");
  }

  Arith::Arith(DOMElement* elt) : m_elt(elt), m_tag(-1), m_evaluated(false) {
    std::string tagNameStr;
    int i = 0;
    bool notFound = true;
    if (typeNames == 0) init();
    
    //    m_elt = elt;
    //    m_evaluated = false;
    //    m_tag = -1;
    if (elt) {

      tagNameStr = xmlBase::Dom::getTagName(elt);

      while ((notFound) && (i < ETAG_n)) {
        ptrString p = typeNames[i];
        if (*p == tagNameStr) {
          m_tag = i;
          notFound = false;
        }
        i++;
      }
    }
    if (notFound) {
      std::string msg = std::string("From Arith::Arith unrecognized tag ") +
        tagNameStr;
      throw BadTag(msg);
    }
  }

  double Arith::evaluate() {
    if (!m_evaluated) { // evaluate
      DOMElement* curElt;
      switch(m_tag) {
      case ETAG_const: {
        if (xmlBase::Dom::hasAttribute(m_elt, valString->c_str())) {
          try {
            m_number = 
              xmlBase::Dom::getDoubleAttribute(m_elt, std::string(*valString));
          }
          catch(xmlBase::WrongAttributeType ex) {
            std::cerr << std::endl << ex.getMsg() << std::endl;
            throw ex;
          }
          m_number *= getUnitsScale(m_elt);
        }
        else  { // must have a single operator child or refer child
          // either of which can be evaluated, optionally preceded
          // by a <notes> child, which we ignore
          curElt = xmlBase::Dom::getFirstChildElement(m_elt);

          if (*notesString == xmlBase::Dom::getTagName(curElt)) {
            curElt = xmlBase::Dom::getSiblingElement(curElt);
          }
          m_number = Arith(curElt).evaluate();
        }
        break;
      }
      case ETAG_refer: {
        // Find the element pointed to and evaluate it
        std::string ref = xmlBase::Dom::getAttribute(m_elt, *refToString);
        DOMDocument *doc = m_elt->getOwnerDocument();
        curElt = xmlBase::Dom::getElementById(doc, ref);

        if (std::string("prim") == xmlBase::Dom::getTagName(curElt) ) {

          try {
            m_number = 
              xmlBase::Dom::getDoubleAttribute(curElt, *valString );
          }
          catch(xmlBase::WrongAttributeType ex) {
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
        curElt = xmlBase::Dom::getFirstChildElement(m_elt);
        m_number = -(Arith(curElt).evaluate());
        break;
      }
      case ETAG_half: {
        // Have a single child
        curElt = xmlBase::Dom::getFirstChildElement(m_elt);
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
        std::string msg = 
          std::string("From Arith::evaluate  Unknown element type");
        throw BadTag(msg);
      }
      }
      m_evaluated = true;
    }
    return m_number;
  }

  double Arith::add() {
    double  ans = 0.0;
    DOMElement* elt = xmlBase::Dom::getFirstChildElement(m_elt);
    ans += Arith(elt).evaluate();
    elt = xmlBase::Dom::getSiblingElement(elt);
    while (elt != 0) {
      ans += Arith(elt).evaluate();
      elt = xmlBase::Dom::getSiblingElement(elt);
    }
    return ans;
  }

  double Arith::mul() {
    double  ans = 1.0;
    DOMElement* elt = xmlBase::Dom::getFirstChildElement(m_elt);
    ans *= Arith(elt).evaluate();
    elt = xmlBase::Dom::getSiblingElement(elt);
    while (elt != 0) {
      ans *= Arith(elt).evaluate();
      elt = xmlBase::Dom::getSiblingElement(elt);
    }
    return ans;
  }

  double Arith::myMax() {
    DOMElement* elt = xmlBase::Dom::getFirstChildElement(m_elt);
    double  ans = Arith(elt).evaluate();
    elt = xmlBase::Dom::getSiblingElement(elt);
    while (elt != 0) {
      double newAns = Arith(elt).evaluate();
      if (newAns > ans) ans = newAns;
      elt = xmlBase::Dom::getSiblingElement(elt);
    }
    return ans;
  }
  double Arith::minus() {
    double  ans;
    DOMElement* elt = xmlBase::Dom::getFirstChildElement(m_elt);
    ans = Arith(elt).evaluate();
    elt = xmlBase::Dom::getSiblingElement(elt);
    ans -= Arith(elt).evaluate();
    return ans;
  }

  double Arith::quo() {
    double  ans;
    double  divisor;
    DOMElement* elt = xmlBase::Dom::getFirstChildElement(m_elt);
    ans = Arith(elt).evaluate();
    elt = xmlBase::Dom::getSiblingElement(elt);
    divisor = Arith(elt).evaluate();  // check for 0?
    ans /= divisor;
    return ans;
  }

  // Store value back in DOM (evaluating first if necessary)
  void Arith::saveValue() {
    // Only makes sense for a const
    if (m_tag != ETAG_const) return;

    if (!m_evaluated) evaluate();

    // Answer will always be returned in mm for length-type
    // constants, so, just in case, this attribute should be
    // set to "mm".  If we're not a length-type constant,
    // setting this attribute is harmless.

    xmlBase::Dom::addAttribute(m_elt, *valString, m_number);
    xmlBase::Dom::addAttribute(m_elt, std::string("unitLength"), 
                           std::string("mm"));
  }

  double Arith::getUnitsScale(const DOMElement* elt) {
    if (xmlBase::Dom::getAttribute(elt, "uType") != *lengthString) {
      return 1.0;
    }

    
    std::string unitLength = xmlBase::Dom::getAttribute(elt, "unitLength");
    if (unitLength == *cmString) return 10.0;
    else if (unitLength == *mString) return 1000.0;
    else return 1.0;
  }


}

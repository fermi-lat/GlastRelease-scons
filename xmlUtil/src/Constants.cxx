// $Header$

#include <string>
#include <xercesc/dom/DOMString.hpp>
#include <xercesc/dom/DOM_Node.hpp>
#include <xercesc/dom/DOM_NodeList.hpp>
#include <xercesc/dom/DOM_TreeWalker.hpp>
#include "xmlUtil/Constants.h"
#include "xmlUtil/Arith.h"
#include "xml/Dom.h"
#include <vector>

namespace {
  // Put a local utility in the unnamed namespace

/* For each <prim>
    if it's of type "length"
     if units are not "mm" then convert to "mm"

   for each <primEnergy>
      Make a <prim> node with same name, value fields to replace it.
      Convert GeV to MeV

   For diagnostics, digits in return value are calculated as follows:
           1's place (1)  -- was a <prim>
              1000 place  -- was a <primEnergy>
       10000 place  -- required conversion
      100000 place  -- unknown unit; couldn't convert
*/

  int normPrim(DOM_Element elt, std::vector<DOM_Element>& save) {
    using xml::Dom;

    int ret = 0;
    if (Dom::checkTagName(elt, "prim") ) {
    
      ret += 1;
      bool isInt = false;

      // Used to check if we're supposed to be an int and  If so, coerce
      // m_number to be nearby int in case of round-off error.  No
      // longer do this.   
      std::string typeAtt = Dom::getAttribute(elt, "type");
      if (isInt = (typeAtt == std::string("int"))) {
        // New style:  it had better be an int!  If not, it's an error
        try {
          int intVal = Dom::getIntAttribute(elt, "value");
        }
        catch (xml::DomException ex) {
          std::cerr << "from xml::Constants::normPrim " << ex.getMsg() 
                    << std::endl;
          throw ex;
        }
      }

      // Found a prim, but we don't have to convert since it's not a length
      if (!(elt.getAttribute(DOMString("uType"))).equals("length")) return ret;
      if (isInt) {   // mistake.  Don't use type=int for length

      }

      std::string units = Dom::getAttribute(elt, "unitLength");
      if (units == std::string("mm")) return ret;
    
      double scale, value;
      if (units == std::string("cm"))  {
        ret += 10000;
        scale = 10;
      }

      else if (units == std::string("m")) {
        ret += 10000;
        scale = 1000;
      }
      else return (100000 + ret);          // should never happen
    
      try {
        value = Dom::getDoubleAttribute(elt, "value");
      }
      catch (xml::DomException ex) {
        std::cerr << "from xml::Constants::normPrim " << ex.getMsg() 
                  << std::endl;
        throw ex;
      }

      value *= scale;
      Dom::addAttribute(elt, std::string("value"), value);
      Dom::addAttribute(elt, std::string("unitLength"), std::string("mm"));
    
      //      return ret;
    }
    else if (Dom::checkTagName(elt, "primEnergy") ) {
      save.push_back(elt);
      ret += 1000;
    }
    return ret;
  }                 // end of utility normPrim

  void handleEnergies(std::vector<DomElement> saved) {
    // Make a new <prim> node to replace it; meanwhile convert
    // GeV to MeV if necessary
    using xml::Dom;

    std::vector<DomElement>::const_iterator eltIt;
    for (eltIt = saved.begin(); eltIt != saved.end(); eltIt++) {
      DomElement elt = *eltIt;
      DOM_Document doc = elt.getOwnerDocument();
      DomElement  prim = doc.createElement("prim");
      Dom::addAttribute(prim, "name", Dom::getAttribute(elt, "name"));
      Dom::addAttribute(prim, "type", "double");
      Dom::addAttribute(prim, "uType", "energy");
      Dom::addAttribute(prim, "unitEnergy", "MeV");

      double energy;
      try {
        energy = Dom::getDoubleAttribute(elt, "value");
      }
      catch (xml::DomException ex) {
        std::cerr << "From Constants;handleEnergies " << ex.getMsg() 
                  << std::endl;
        throw ex;
      }

      if (std::string("GeV") == Dom::getAttribute(elt, "units") ) {
        energy *= 1000;
      }

      Dom::addAttribute(prim, std::string("value"), energy);
      elt.normalize();
      DomNode textChild = elt.getFirstChild();
      prim.appendChild(textChild);
      //      elt.removeChild(textChild);
      DomNode parent = elt.getParentNode();
      DomNode oldChild =  parent.replaceChild(prim, elt);
    }
  }
}   // end of anonymous namespace

namespace xmlUtil {

  Constants::Constants(DOM_Document doc) : m_doc(doc) {

    // There is at most 1 <constants> element
    m_constants = 
      xml::Dom::findFirstChildByName(doc.getDocumentElement(), "constants");
  }

  void Constants::normalizePrimary()
  {
    using xml::Dom;

    if (m_constants == DomElement()) return;

    DomElement primary = Dom::findFirstChildByName(m_constants, "primary");

    DomNode curNode = primary;

    if (std::string("true") == Dom::getAttribute(primary, "normalized"))
      return;

    // Make a DOMTreeWalker which returns only elements
    unsigned long whatToShow = 1 << (DomNode::ELEMENT_NODE -1);
    
    DOM_TreeWalker walker = 
      //      m_doc.createTreeWalker(curNode, whatToShow, 0, 0);
      m_doc.createTreeWalker(primary, whatToShow, 0, 0);
    
    // Diagnostic -- number of <prim>s seen
    int count = 0;
    std::vector<DomElement> energyElts;
    while (curNode != DomNode() ) {
      count += normPrim(static_cast<DomElement &> (curNode), energyElts);
      curNode = walker.nextNode();
    }
    handleEnergies(energyElts);

    Dom::addAttribute(primary, "normalized", "true");
  }
 
  // Remove children of evaluated constants
  void Constants::pruneConstants(bool keepNotes) { 
    using xml::Dom;

    if (m_constants == DomNode()) return;

    DomElement curConst;
    
    // For each derCategory, prune each const child
    std::vector<DomElement> cats;

    Dom::getDescendantsByTagName(m_constants, "derCategory", cats);
    int nCats = cats.size();
    int iCat;
    
    for (iCat = 0; iCat < nCats; iCat++) {
      DomElement& curCat = cats[iCat];
      curConst = Dom::findFirstChildByName(curCat, "const" );
      
      while (curConst != DomElement()) {
        if (keepNotes) { // tread more carefully
          DomElement child = Dom::getFirstChildElement(curConst);
          if (child != DomElement() ) {

            if (Dom::checkTagName(child, "notes") ) {
              child = Dom::getSiblingElement(child);
              if (child != DomElement()) {
                Dom::prune(child);
                (child.getParentNode()).removeChild(child);
              }
            }
            else {
              Dom::prune(curConst);  // remove all its children
            }
          }
        }
        else {  // remove all its children
          Dom::prune(curConst);  
        }
        curConst = Dom::getSiblingElement(curConst);
      }
    }
  }

  /*! For each derCategory, evaluate each const child, writing
      value back to "value" attribute in the DOM */
  void Constants::evalConstants() {
    using xml::Dom;

    if (m_constants == DomNode()) return;
    DomElement derived = 
      Dom::findFirstChildByName(m_constants, "derived"); 
    if (derived == DomElement()) return; 
    if (std::string("true") == Dom::getAttribute(derived, "evaluated")) return;

    DomElement curConst;

    std::vector<DomElement> cats;
    
    Dom::getDescendantsByTagName(m_doc.getDocumentElement(), "derCategory",
                              cats);
    int nCats = cats.size();
    int iCat;
    
    for (iCat = 0; iCat < nCats; iCat++) {
      DomElement& curCat = cats[iCat];
      curConst = Dom::findFirstChildByName(curCat, "const" );
      
      while (curConst != DomElement()) {
        try {
          Arith curArith(curConst);
          double evalValue = curArith.evaluate();
          curArith.saveValue();
        }
        catch (XmlUtilException ex) {
          std::cerr << ex.getMsg();
          throw(ex);
        }
        curConst = xml::Dom::getSiblingElement(curConst);
      }
    }
    Dom::addAttribute(derived, "evaluated", "true");
  }
}


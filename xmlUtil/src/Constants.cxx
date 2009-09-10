// $Header$
#include <iostream>
#include <string>

#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMDocument.hpp>

#include <xercesc/dom/DOMTreeWalker.hpp>
#include <xercesc/util/XMLString.hpp>
#include "xmlUtil/Constants.h"
#include "xmlUtil/Arith.h"
#include "xmlBase/Dom.h"

namespace {
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMNode;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMNodeFilter;
  using XERCES_CPP_NAMESPACE_QUALIFIER XMLString;

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

  int normPrim(DOMElement* elt, std::vector<DOMElement*>& save) {
    using xmlBase::Dom;

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
          /*int intVal = */
          Dom::getIntAttribute(elt, "value");
        }
        catch (xmlBase::DomException ex) {
          std::cerr << "from xmlBase::Constants::normPrim " << ex.getMsg() 
                    << std::endl;
          throw ex;
        }
      }

      // Found a prim, but we don't have to convert since it's not a length
      // if (!(elt.getAttribute(DOMString("uType"))).equals("length")) return ret;
      if (!(xmlBase::Dom::getAttribute(elt, "uType") == "length" ) ) return ret;

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
      catch (xmlBase::DomException ex) {
        std::cerr << "from xmlBase::Constants::normPrim " << ex.getMsg() 
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

}   // end of anonymous namespace

namespace xmlUtil {
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMNode;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMTreeWalker;

  Constants::Constants(DOMDocument* doc) : m_doc(doc), m_constants(0),
                                           m_walker(0) {
    // There is at most 1 <constants> element
    m_constants = 
      xmlBase::Dom::findFirstChildByName(doc->getDocumentElement(), "constants");
    if (m_constants) {
      m_walker = doc->createTreeWalker(m_constants, 
                                       DOMNodeFilter::SHOW_ELEMENT,
                                       0, true);
    }
  }

  // Deleting the TreeWalker causes free(  ) to complain.  Could
  // be that deletion is handled by the parent document
  Constants::~Constants() { /* delete m_walker; */} 

  void Constants::handleEnergies(std::vector<DOMElement*> saved) {
    // Make a new <prim> node to replace it; meanwhile convert
    // GeV to MeV if necessary
    using xmlBase::Dom;

    std::vector<DOMElement*>::const_iterator eltIt;
    for (eltIt = saved.begin(); eltIt != saved.end(); eltIt++) {
      DOMElement* elt = *eltIt;
      DOMDocument* doc = elt->getOwnerDocument();
      XMLCh* xmlchPrim = XMLString::transcode("prim");
      DOMElement*  prim = doc->createElement(xmlchPrim);
      XMLString::release(&xmlchPrim);
      Dom::addAttribute(prim, "name", Dom::getAttribute(elt, "name"));
      Dom::addAttribute(prim, "type", "double");
      Dom::addAttribute(prim, "uType", "energy");
      Dom::addAttribute(prim, "unitEnergy", "MeV");

      double energy;
      try {
        energy = Dom::getDoubleAttribute(elt, "value");
      }
      catch (xmlBase::DomException ex) {
        std::cerr << "From Constants;handleEnergies " << ex.getMsg() 
                  << std::endl;
        throw ex;
      }

      if (std::string("GeV") == Dom::getAttribute(elt, "units") ) {
        energy *= 1000;
      }

      Dom::addAttribute(prim, std::string("value"), energy);
      elt->normalize();
      DOMNode* textChild = elt->getFirstChild();
      prim->appendChild(textChild);
      //      elt.removeChild(textChild);
      //      DOMNode* parent = elt->getParentNode();
      m_walker->setCurrentNode(elt);
      DOMNode* parent = m_walker->parentNode();
      /*DOMNode* oldChild =  */
      parent->replaceChild(prim, elt);
    }
  }

  void Constants::normalizePrimary()
  {
    using xmlBase::Dom;

    if (m_constants == 0) return;

    DOMElement* primary = Dom::findFirstChildByName(m_constants, "primary");

    DOMNode* curNode = primary;

    if (std::string("true") == Dom::getAttribute(primary, "normalized"))
      return;

    // Make a DOMTreeWalker which returns only elements
    /*  Now we already have one
    unsigned long whatToShow = 1 << (DOMNode::ELEMENT_NODE -1);
    
    DOMTreeWalker* walker = 
      //      m_doc.createTreeWalker(curNode, whatToShow, 0, 0);
      //     last argument says to expand entity references
      m_doc->createTreeWalker(primary, whatToShow, 0, true);
    */
    m_walker->setCurrentNode(curNode);
    // Diagnostic -- number of <prim>s seen
    int count = 0;
    std::vector<DOMElement*> energyElts;
    while (curNode != 0 ) {
      count += normPrim(static_cast<DOMElement*> (curNode), energyElts);
      curNode = m_walker->nextNode();
    }
    handleEnergies(energyElts);

    Dom::addAttribute(primary, "normalized", "true");
  }
 
  // Remove children of evaluated constants
  void Constants::pruneConstants(bool keepNotes) { 
    using xmlBase::Dom;

    if (m_constants == 0) return;

    DOMElement* curConst;
    
    // For each derCategory, prune each const child
    std::vector<DOMElement*> cats;

    Dom::getDescendantsByTagName(m_constants, "derCategory", cats);
    int nCats = cats.size();
    int iCat;
    
    for (iCat = 0; iCat < nCats; iCat++) {
      DOMElement* curCat = cats[iCat];
      curConst = Dom::findFirstChildByName(curCat, "const" );
      
      while (curConst != 0) {
        if (keepNotes) { // tread more carefully
          DOMElement* child = Dom::getFirstChildElement(curConst);
          if (child != 0 ) {

            if (Dom::checkTagName(child, "notes") ) {
              child = Dom::getSiblingElement(child);
              if (child != 0) {
                Dom::prune(child);
                (child->getParentNode())->removeChild(child);
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
    using xmlBase::Dom;

    if (m_constants == 0) return;
    DOMElement* derived = 
      Dom::findFirstChildByName(m_constants, "derived"); 
    if (derived == 0) return; 
    if (std::string("true") == Dom::getAttribute(derived, "evaluated")) return;

    DOMElement* curConst;

    std::vector<DOMElement*> cats;
    
    Dom::getDescendantsByTagName(m_doc->getDocumentElement(), "derCategory",
                              cats);
    int nCats = cats.size();
    int iCat;
    
    for (iCat = 0; iCat < nCats; iCat++) {
      DOMElement* curCat = cats[iCat];
      curConst = Dom::findFirstChildByName(curCat, "const" );
      while (curConst != 0) {
        try {
          Arith curArith(curConst);
          /* double evalValue = */
          curArith.evaluate();

          /*
                for debug
          std::string curConstName = Dom::getAttribute(curConst, "name");
          DOMElement* eltById = Dom::getElementById(m_doc, curConstName);
          XMLCh* xmlchConstName = XMLString::transcode(curConstName.c_str());
          DOMElement* XercesEltById = m_doc->getElementById(xmlchConstName);
          std::cout << "curConst addr: " << curConst << " my ById: " 
                    << eltById <<  " Xerces byId: " << xmlchConstName << std::endl;
          XMLString::release(&xmlchConstName);
          */
          curArith.saveValue();
          /*

            TEMPORARILY eliminate "save" to see if this makes a difference

          debug code
          if (curConstName == "bigXZDim") {
            curArith.saveValue();
            std::cout << curConstName << " evaluated and saved " << std::endl;
            std::string bigXZDimValue = Dom::getAttribute(curConst, "value");
            std::cout << "..saved value is " << bigXZDimValue << std::endl;
          }
          else curArith.saveValue();
          */
        }
        catch (XmlUtilException ex) {
          std::cerr << ex.getMsg();
          throw(ex);
        }
        curConst = xmlBase::Dom::getSiblingElement(curConst);
      }
    }
    Dom::addAttribute(derived, "evaluated", "true");
  }
}


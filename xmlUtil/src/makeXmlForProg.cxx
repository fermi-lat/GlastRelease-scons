// $Header$
/*! \file Standalone program to transform source xml file into a
    preprocessed version suitable for most clients programs (such
    as Simulation and Reconstruction).     Clients needing to
    extract documentation from the xml file will go a different
    route (see makeXmlForDoc.cxx.)

    In particular this program will
      - add a <source> element to the output, indicating how and
        from what the output was generated
      - convert any <prim> values not already expressed in standard
        units (mm) to the standard
      - substitute values for all references constants, evaluating
        such constants as needed
      - delete any <derCategory> without the "save" attribute set
        to true and all its children.  For any <derCategory> which
        is to be saved, delete all element content from each
        <const> (except for notes?).
 */

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "xmlUtil/Source.h"
#include "xmlUtil/Constants.h"
#include "xmlUtil/Arith.h"
#include "xmlUtil/Substitute.h"
#include <xercesc/dom/DOM_Element.hpp>
#include <xercesc/dom/DOM_NodeList.hpp>
#include <xercesc/dom/DOM_NamedNodeMap.hpp>
#include <xercesc/dom/DOM_DocumentType.hpp>

#include <string>
#include <iostream>
#include <fstream>

std::ostream *openOut(char * outfile);
void outProlog(const DOM_DocumentType& doctype, std::ostream& out);

const std::string myId("$Id$");

// Can't literally put in the string we want or CVS will mess it up.
// Instead make a copy of this template, replacing the # with $
const std::string idTemplate("#Id: not committed $");

/*!
    Main program for the eval application.
    \param arg1 is the input xml file
    \param arg2 is specification for the output file, or "-" to output
    to std::cout
*/
int main(int argc, char* argv[]) {
  std::ostream *out;
  if (argc < 3) {  // instructions
    std::cout << "Required first argument is xml file to be parsed" 
              << std::endl;
    std::cout << "Required second argument is output file (- for stdout)"
              << std::endl;
    exit(0);
  }

  xml::XmlParser* parser = new xml::XmlParser();
  DOM_Document doc = parser->parse(argv[1], "gdd");

  if (doc == 0) {
    std::cout << "Document failed to parse correctly" << std::endl;
    exit(0);
  }

  // Else successful.   Open output
  out = openOut(argv[2]);

  DOM_Element docElt = doc.getDocumentElement();
  DOM_DocumentType  doctype = doc.getDoctype(); //  check for gdd?? 
  xmlUtil::Constants   *constants = new xmlUtil::Constants(doc);

  // prim processing
  constants->normalizePrimary();

  // Evaluate all constants.  Most would be evaluated even without
  // this step because of the use of Substitute, but if there happens
  // to be a constant not actually used in the <section> part of
  // the file, Substitute would not cause that constant to be evaluated.
  constants->evalConstants();


  xmlUtil::Substitute* sub = new xmlUtil::Substitute(doc);

  std::vector<DOM_Element> dicts;
  xml::Dom::getDescendantsByTagName(docElt, "idDict", dicts);

  //  There is only one dictionary
  DOM_Element& dictElt = dicts[0];
  int nDictSub = sub->execute(dictElt);
    std::cout << "#elements substituted for id dictionary:  " << nDictSub
              << std::endl;

  // Search for all dimension and position references, substitute
  // value of referenced element.
  std::vector<DOM_Element> sections;
  xml::Dom::getDescendantsByTagName(docElt, "section", sections);

  unsigned int nSec = sections.size();
  unsigned int iSec;
  for (iSec = 0; iSec < nSec; iSec++) {
    int nSub;
    DOM_Element& secElt = sections[iSec];
    nSub = sub->execute(secElt);
    std::cout << "#elements substituted for in this section:  " << nSub
              << std::endl;
  }
  
  // Get rid of any <derCategory> without the save attribute set to "true".
  // For the rest, prune element content from all <const> children;
  // it's no longer needed.
  DOM_Element derived = xml::Dom::findFirstChildByName(docElt, "constants");
  if (derived != DOM_Element() ) {
    derived = xml::Dom::findFirstChildByName(derived, "derived");
    
    if (derived != DOM_Element()) {

      std::vector<DOM_Element> cats;
      xml::Dom::getDescendantsByTagName(docElt, "derCategory", cats);

      unsigned nCat = cats.size();

      //      DOM_NodeList cats = derived.getElementsByTagName("derCategory");
      //      int nCat = cats.getLength();
      unsigned int iCat;
      std::vector<DomElement> toRemove;
      for (iCat=0; iCat < nCat; iCat++) {
        DOM_Element elt = cats[iCat];
        if (std::string("false") ==  xml::Dom::getAttribute(elt, "save") ) {
          //remove
        
          xml::Dom::prune(elt);
          toRemove.push_back(elt);
        }
      }
      if (unsigned nRemove = toRemove.size() > 0) {
        // all elements to be removed have the same parent
        DOM_Node parent = (toRemove[0]).getParentNode();
        for (unsigned iVictim= 0; iVictim < nRemove; iVictim++) {
          parent.removeChild(toRemove[iVictim]);
        }
      }
      // Now just get rid of all content for <const> children 
      // in remaining derived categories
      constants->pruneConstants();
    }
  }

  // Add a <source> child to the outer gdd element
  xmlUtil::Source * source = 
    new xmlUtil::Source(doc, "xmlUtil/v1/src/forProg.exe", myId.c_str());
  source->add();
  
  // Output the xml declaration and all the text in the DOCTYPE (see DOMPrint)
  outProlog(doctype, *out);

  // If have gdd element with CVSid attribute, null it out.  Don't have
  // a real CVS id until the file has been committed
  if (docElt.getAttribute("CVSid") != DOMString() ) {
    std::string noId(idTemplate);
    noId.replace(0, 1, "$");
    xml::Dom::addAttribute(docElt, "CVSid", noId);
    //    docElt.setAttribute("CVSid", noId.c_str());
  }

  // Finally output the elements
  // May want option to exclude comments here
  xml::Dom::prettyPrintElement(docElt, *out, "");

  delete parser;
  return(0);
}

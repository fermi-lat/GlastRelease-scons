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

#include "xmlBase/XmlParser.h"
#include "xmlBase/Dom.h"
#include "xmlUtil/Source.h"
#include "xmlUtil/Constants.h"
#include "xmlUtil/Arith.h"
#include "xmlUtil/Substitute.h"
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMDocumentType.hpp>

#include <string>
#include <iostream>
#include <fstream>

std::ostream *openOut(char * outfile);
void outProlog(const XERCES_CPP_NAMESPACE_QUALIFIER DOMDocumentType* doctype, 
               std::ostream& out);

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
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocumentType;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMNode;

  std::ostream *out;
  if (argc < 3) {  // instructions
    std::cout << "Required first argument is xml file to be parsed" 
              << std::endl;
    std::cout << "Required second argument is output file (- for stdout)"
              << std::endl;
    exit(0);
  }

  xmlBase::XmlParser* parser = new xmlBase::XmlParser();
  DOMDocument* doc = parser->parse(argv[1], "gdd");

  if (doc == 0) {
    std::cout << "Document failed to parse correctly" << std::endl;
    exit(0);
  }

  // Else successful.   Open output
  out = openOut(argv[2]);

  DOMElement* docElt = doc->getDocumentElement();
  DOMDocumentType*  doctype = doc->getDoctype(); //  check for gdd?? 
  xmlUtil::Constants   *constants = new xmlUtil::Constants(doc);

  // prim processing
  constants->normalizePrimary();

  // Evaluate all constants.  Most would be evaluated even without
  // this step because of the use of Substitute, but if there happens
  // to be a constant not actually used in the <section> part of
  // the file, Substitute would not cause that constant to be evaluated.
  constants->evalConstants();


  xmlUtil::Substitute* sub = new xmlUtil::Substitute(doc);

  std::vector<DOMElement*> dicts;
  xmlBase::Dom::getDescendantsByTagName(docElt, "idDict", dicts);

  //  There is only one dictionary
  DOMElement* dictElt = dicts[0];
  int nDictSub = sub->execute(dictElt);
    std::cout << "#elements substituted for id dictionary:  " << nDictSub
              << std::endl;

  // Search for all dimension and position references, substitute
  // value of referenced element.
  std::vector<DOMElement*> sections;
  xmlBase::Dom::getDescendantsByTagName(docElt, "section", sections);

  unsigned int nSec = sections.size();
  unsigned int iSec;
  for (iSec = 0; iSec < nSec; iSec++) {
    int nSub;
    DOMElement* secElt = sections[iSec];
    nSub = sub->execute(secElt);
    std::cout << "#elements substituted for in this section:  " << nSub
              << std::endl;
  }
  
  // Get rid of any <derCategory> without the save attribute set to "true".
  // For the rest, prune element content from all <const> children;
  // it's no longer needed.
  DOMElement* derived = xmlBase::Dom::findFirstChildByName(docElt, "constants");
  if (derived != 0 ) {
    derived = xmlBase::Dom::findFirstChildByName(derived, "derived");
    
    if (derived != 0) {

      std::vector<DOMElement*> cats;
      xmlBase::Dom::getDescendantsByTagName(docElt, "derCategory", cats);

      unsigned nCat = cats.size();

      unsigned int iCat;
      std::vector<DOMElement*> toRemove;
      for (iCat=0; iCat < nCat; iCat++) {
        DOMElement* elt = cats[iCat];
        if (std::string("false") ==  xmlBase::Dom::getAttribute(elt, "save") ) {
          //remove
        
          xmlBase::Dom::prune(elt);
          toRemove.push_back(elt);
        }
      }
      if (unsigned nRemove = toRemove.size() > 0) {
        // all elements to be removed have the same parent
        DOMNode* parent = (toRemove[0])->getParentNode();
        for (unsigned iVictim= 0; iVictim < nRemove; iVictim++) {
          parent->removeChild(toRemove[iVictim]);
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
  //  if (docElt->getAttribute("CVSid") != 0 ) {
  if (xmlBase::Dom::hasAttribute(docElt, "CVSid")) {
    std::string noId(idTemplate);
    noId.replace(0, 1, "$");
    xmlBase::Dom::addAttribute(docElt, "CVSid", noId);
    //    docElt.setAttribute("CVSid", noId.c_str());
  }

  // Finally output the elements
  // May want option to exclude comments here
  xmlBase::Dom::prettyPrintElement(docElt, *out, "");

  delete parser;
  return(0);
}

// $Header$
/*! \file Standalone program to transform source xml file into a 
    preprocessed version suitable for documentation.  Typically
    will be transformed further, e.g. to html by an xslt transform.
    xslt is not suitable to do the whole job because it can't handle
    recursive evaluation of the constants.
  
    In particular this program will
      - add a <source> element to the output, indicating how and
        from what the output was generated
      - evaluate all <const>
      - delete all content except <notes> elements from each <const>
      - delete all <section> elements (for now.  Later may think of
        some useful information for doc. purposes to be extracted
        from them).
      - delete all materials 
      - output the resulting xml document.
 */

#include "xmlBase/XmlParser.h"
#include "xmlBase/Dom.h"
#include "xmlUtil/Source.h"
#include "xmlUtil/Arith.h"
#include "xmlUtil/Constants.h"
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
    Main program for the forDoc application.
    \param arg1 is the input xml file
    \param arg2 is specification for the output file, or "-" to output
    to std::cout
*/
int main(int argc, char* argv[]) {

  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocumentType;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

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

  // For each const child of a derCategory
  //   evaluate it, saving value
  //   delete all children except <notes>
  xmlUtil::Constants *constants = new xmlUtil::Constants(doc);
  constants->evalConstants();
  constants->pruneConstants(true);  

  // Delete any id dictionaries

  std::vector<DOMElement*> dicts;
  xmlBase::Dom::getDescendantsByTagName(docElt, "idDict", dicts);

  unsigned nDict = dicts.size();
  for (unsigned iDict = 0; iDict < nDict; iDict++) {
    DOMElement* dictElt = dicts[iDict];

    xmlBase::Dom::prune(dictElt);
    (dictElt->getParentNode())->removeChild(dictElt);
  }

  // Delete all sections
  std::vector<DOMElement*> sections;
  xmlBase::Dom::getDescendantsByTagName(docElt, "section", sections);

  unsigned nSec = sections.size();
  for (unsigned iSec = 0; iSec < nSec; iSec++) {

    DOMElement* secElt = sections[iSec];
    xmlBase::Dom::prune(secElt);
    (secElt->getParentNode())->removeChild(secElt);
  }

  // Delete materials
  std::vector<DOMElement*> materials;
  xmlBase::Dom::getDescendantsByTagName(docElt, "materials", materials);

  unsigned nM = materials.size();
  for (unsigned iM = 0; iM < nM; iM++) {

    DOMElement* matElt = materials[iM];
    xmlBase::Dom::prune(matElt);
    (matElt->getParentNode())->removeChild(matElt);
    //    secNode = toCome;
  }

  // Add a <source> child to the outer gdd element
  xmlUtil::Source * source = 
    new xmlUtil::Source(doc, "xmlUtil/v1/src/forDoc.exe", myId.c_str());
  source->add();
  
  // Output the xml declaration and all the text in the DOCTYPE (see DOMPrint)
  outProlog(doctype, *out);

  // If have gdd element with CVSid attribute, null it out.  Don't have
  // a real CVS id until the file has been committed
  std::string cvsId = xmlBase::Dom::getAttribute(docElt, "CVSid");
  if (cvsId.size() > 0) {
    std::string noId(idTemplate);
    noId.replace(0, 1, "$");
    xmlBase::Dom::addAttribute(docElt, "CVSid", noId);
  }
  // Finally output the elements
  // May want option to exclude comments here
  xmlBase::Dom::prettyPrintElement(docElt, *out, "");

  delete parser;
  return(0);
}

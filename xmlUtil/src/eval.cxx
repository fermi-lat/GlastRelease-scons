// $Header$
/*! \file Standalone program to transform source xml file into a 
    preprocessed version suitable for most clients (documentation-type
    clients will probably stick with the original file). 
  
    In particular this program will
      - add a <source> element to the output, indicating how and
        from what the output was generated
      - convert any <prim> values not already expressed in standard
        units (mm) to the standard
      - substitute values for all references constants, evaluating
        such constants as needed
      - delete the <derived> element and all descendents from the
        output
 */

#include "xmlBase/XmlParser.h"
#include "xmlBase/Dom.h"
#include "xmlUtil/Arith.h"
#include "xmlUtil/Substitute.h"
#include "xmlUtil/Constants.h"
#include "xmlUtil/Source.h"
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMDocumentType.hpp>

#include <string>
#include <iostream>
#include <fstream>

// openOut and outProlog can be found in src/local/outUtils.cxx
std::ostream *openOut(char * outfile);
void outProlog(const XERCES_CPP_NAMESPACE_QUALIFIER DOMDocumentType* doctype, 
               std::ostream& out);

const std::string myId("$Id$");

/*!
    Main program for the eval application.
    \param arg1 is the input xml file
    \param arg2 is specification for the output file, or "-" to output
    to std::cout
*/
int main(int argc, char* argv[]) {
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocumentType;

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

  // prim processing
  xmlUtil::Constants *constants = new xmlUtil::Constants(doc);
  constants->normalizePrimary();

  // No longer need the evaluate step since substitution will
  // evaluate any constant which has not yet been evaluated.

  // Search for all dimension and position references, substitute
  // value of referenced element.

  std::vector<DOMElement*> sections;
  xmlBase::Dom::getChildrenByTagName(docElt, "section", sections);
  int nSec = sections.size();
  int iSec;

  xmlUtil::Substitute* sub = new xmlUtil::Substitute(doc);

  for (iSec = 0; iSec < nSec; iSec++) {
    int nSub;
    DOMElement* secElt = sections[iSec];
    nSub = sub->execute(secElt);
    std::cout << "#elements substituted for in this section:  " << nSub
              << std::endl;
  }

  // In case there were no sections to substitute, do eval here
  constants->evalConstants();


  // Add a <source> child to the outer gdd element
  xmlUtil::Source *source = 
    new xmlUtil::Source(doc, "xmlUtil/v1/src/eval.exe", "$Id$");
  source->add();
  
  // Output the xml declaration and all the text in the DOCTYPE (see DOMPrint)
  outProlog(doctype, *out);

  // Finally output the elements
  // May want option to exclude comments here
  xmlBase::Dom::prettyPrintElement(docElt, *out, "");

  delete parser;
  return(0);
}








// $Header$
/*! \file Standalone program to test Arith class

 */

#include "xmlBase/XmlParser.h"
#include "xmlBase/Dom.h"
#include "xmlUtil/Arith.h"
#include "xmlUtil/Substitute.h"
#include "xmlUtil/Constants.h"
#include "xmlUtil/Source.h"
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMDocumentType.hpp>
#include <xercesc/util/XMLString.hpp>



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
    \param arg2 is specification for the output file, or "-" to output -- so 
    to std::cout
*/
int main(int argc, char* argv[]) {
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMNode;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMNodeList;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocumentType;
  using XERCES_CPP_NAMESPACE_QUALIFIER XMLString;

  // std::ostream *out;
  if (argc < 2) {  // instructions
    std::cout << "Required first argument is xml file to be parsed" 
              << std::endl;
    //    std::cout << "Required second argument is output file (- for stdout)"
    //              << std::endl;
    exit(0);
  }

  xmlBase::XmlParser* parser = new xmlBase::XmlParser();
  //  DOMDocument* doc = parser->parse(argv[1], "gdd");
  DOMDocument* doc = parser->parse(argv[1], "");

  if (doc == 0) {
    std::cout << "Document failed to parse correctly" << std::endl;
    exit(0);
  }

  // Else successful.   Open output
  //  out = openOut(argv[2]);

  //  DOMElement* docElt = doc->getDocumentElement();
  //  DOMDocumentType*  doctype = doc->getDoctype(); //  check for gdd?? 

  // Just try something very simple, using minimal code outside Xerces,
  //  to see if getElementById is screwed up 

  // Lots of elements with tagname = "const" have ID attribute called "name"
  XMLCh* xmlchConst = XMLString::transcode("const");
  XMLCh* xmlchName = XMLString::transcode("name");
  DOMNodeList* constElts = doc->getElementsByTagName(xmlchConst);

  // For each one, write out address. Get value of name attribute, if any.
  // If has a name, find it via getElementById
  // and write out that address.  Should match
  unsigned int nElt = constElts->getLength();
  for (unsigned int iElt = 0; iElt < nElt; iElt++) {
    DOMNode* item = constElts->item(iElt);
    DOMElement* itemElt = static_cast<DOMElement *>(item);
    std::cout << std::endl << "Const elt " << iElt << " Address as node:  " 
              << item << " and as element: " << itemElt << std::endl;
    const XMLCh* xmlchNamevalue = itemElt->getAttribute(xmlchName);
    if (XMLString::stringLen(xmlchNamevalue) > 0 ) {
      char* namevalue = XMLString::transcode(xmlchNamevalue);
      std::cout << "element has name " << namevalue << std::endl;
      DOMElement* byIdElt = doc->getElementById(xmlchNamevalue);
      std::cout << "Address from getElementById: " << byIdElt << std::endl;
      XMLString::release(&namevalue);
    }
    
  }


  /*
  // prim processing
  DOMElement* bigXZDimElt = xmlBase::Dom::getElementById(doc, "bigXZDim");
  std::cout << "address of const bigXZDim " <<  bigXZDimElt << std::endl;

  std::string XZval = xmlBase::Dom::getAttribute(bigXZDimElt, "value");
  std::cout << "..and its value attr (not expected at this stage " 
            << XZval << std::endl;

  xmlUtil::Constants *constants = new xmlUtil::Constants(doc);
  constants->evalConstants();
  // constants->normalizePrimary();
  bigXZDimElt = xmlBase::Dom::getElementById(doc, "bigXZDim");
  std::cout << "address after eval of const bigXZDim " 
            <<  bigXZDimElt << std::endl;

  XZval = xmlBase::Dom::getAttribute(bigXZDimElt, "value");
  std::cout << "..and its value attr after having done evalConstants "
            << XZval << std::endl;


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
  */
  // In case there were no sections to substitute, do eval here
  //  constants->evalConstants();

  /*
  // Add a <source> child to the outer gdd element
  xmlUtil::Source *source = 
    new xmlUtil::Source(doc, "xmlUtil/v1/src/eval.exe", "$Id$");
  source->add();
  
  // Output the xml declaration and all the text in the DOCTYPE (see DOMPrint)
  outProlog(doctype, *out);

  // Finally output the elements
  // May want option to exclude comments here
  xmlBase::Dom::prettyPrintElement(docElt, *out, "");
  */
  delete parser;
  return(0);
}








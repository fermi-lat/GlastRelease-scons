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

/*!
    Main program for testId application
    \param arg1 is the input xml file
*/
int main(int argc, char* argv[]) {
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMNode;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMNodeList;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocumentType;
  using XERCES_CPP_NAMESPACE_QUALIFIER XMLString;

  //  std::ostream *out;
  if (argc < 2) {  // instructions
    std::cout << "Required first argument is xml file to be parsed" 
              << std::endl;
    exit(0);
  }

  xmlBase::XmlParser* parser = new xmlBase::XmlParser();
  DOMDocument* doc = parser->parse(argv[1], "gdd");

  if (doc == 0) {
    std::cout << "Document failed to parse correctly" << std::endl;
    exit(0);
  }

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
  delete parser;
  return(0);
}








// $Header$

/*! \file Stand-alone test program for id dictionary code */

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "xmlUtil/Substitute.h"
#include "xmlUtil/id/DictValidVisitor.h"
#include "xmlUtil/id/IdDict.h"

#include "dom/DOM_Element.hpp"
#include "dom/DOM_DocumentType.hpp"

#include <iostream>
// #include <fstream>
#include <string>

std::ostream *openOut(char * outfile);

/*!
    Main program for the eval application.
    \param arg1 is the input xml file
*/
int main(int argc, char* argv[]) {
  if (argc < 2) {  // instructions
    std::cout << 
      "Required first argument is xml file (id dictionary) to be parsed" 
              << std::endl;
    exit(0);
  }

  xml::XmlParser* parser = new xml::XmlParser();
  DOM_Document doc = parser->parse(argv[1]);
  DOM_Element  docElt = doc.getDocumentElement();

  if (doc == 0) {
    std::cout << "Document failed to parse correctly" << std::endl;
    exit(0);
  }

  std::cout << "Document was parsed" << std::endl;

  // Evaluate and substitute in case constants were used
  xmlUtil::Substitute* sub = new xmlUtil::Substitute(doc);

  // Look for IdDict element
  DOM_Element dictElt = xml::Dom::findFirstChildByName(docElt, "idDict");

  if (dictElt == DOM_Element()) {
    std::cout << "No identifier dictionary found " << std::endl;
    exit(0);
  }

  std::cout << "Dictionary element found" << std::endl;

  int nSub = sub->execute(dictElt);

  std::cout << " # elements substituted: " << nSub << std::endl;
  
  // Make the dictionary
  xmlUtil::IdDict  *dict = new xmlUtil::IdDict(dictElt);

  std::cout << "Dictionary made" << std::endl;

  // Finally attempt to verify it

  bool valid = dict->isValid();

  std::cout << "Dictionary is " << ( (valid) ? "VALID" : "INVALID" );
  std::cout << std::endl;

  return;
}




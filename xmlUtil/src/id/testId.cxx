// $Header$

/*! \file Stand-alone test program for id dictionary code */

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "xmlUtil/Substitute.h"
#include "xmlUtil/id/DictValidVisitor.h"
#include "xmlUtil/id/IdDict.h"
#include "xmlUtil/id/NamedId.h"

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

  // Test query functions

  xmlUtil::Identifier okId(3);
  xmlUtil::Identifier badId(4);

  okId[0] = 0;
  okId[1] = 1;
  okId[2] = 1;

  badId[0] = 0;
  badId[1] = 0;
  badId[2] = 27;
  badId[3] = 1;

  valid = dict->idOk(okId);
  std::cout << "okId is " <<  ( (valid) ? "VALID" : "INVALID" );
  std::cout << std::endl;

  valid = dict->idOk(badId);
  std::cout << "badId is " <<  ( (valid) ? "VALID" : "INVALID" );
  std::cout << std::endl;

  xmlUtil::NamedId nIdBadVal(5); 
  xmlUtil::NamedId nIdOk(3);
  xmlUtil::NamedId nIdBadName(4);

  nIdBadVal.addField("fTopObjects", 0);
  nIdBadVal.addField("fLATObjects", 0);
  nIdBadVal.addField("fTowerY", 2);
  nIdBadVal.addField("fTowerX", 14);
  nIdBadVal.addField("fTowerObjects", 0);

  nIdOk.addField("fTopObjects", 0);
  nIdOk.addField("fLATObjects", 0);
  nIdOk.addField("fTowerY", 2);

  nIdBadName.addField("fTopObjects", 0);
  nIdBadName.addField("fLATObjects", 0);
  nIdBadName.addField("fTower", 2);
  nIdBadName.addField("fTowerX", 1);

  valid = dict->idOk(nIdBadVal);
  std::cout << "nIdBadVal is " <<  ( (valid) ? "VALID" : "INVALID" );
  std::cout << std::endl;

  valid = dict->idOk(nIdOk);
  std::cout << "nIdOk is " <<  ( (valid) ? "VALID" : "INVALID" );
  std::cout << std::endl;

  valid = dict->idOk(nIdBadName);
  std::cout << "nIdBadNamae is " <<  ( (valid) ? "VALID" : "INVALID" );
  std::cout << std::endl;

  return;
}




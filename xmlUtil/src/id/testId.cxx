// $Header$

/*! \file Stand-alone test program for id dictionary code */

#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "xmlUtil/Substitute.h"
#include "xmlUtil/id/DictValidVisitor.h"
#include "xmlUtil/id/IdDict.h"
#include "xmlUtil/id/NamedId.h"
#include "xmlUtil/id/IdDictMan.h"
#include "xmlUtil/id/IdConverter.h"
#include "xmlUtil/id/IdConversion.h"  // maybe don't need this

#include "dom/DOM_Element.hpp"
#include "dom/DOM_DocumentType.hpp"

#include <iostream>
// #include <fstream>
#include <string>

std::ostream *openOut(char * outfile);

void testQuery(xmlUtil::IdDict* dict);

int testConverter(xmlUtil::IdConverter* converter);

/*!
    Main program for the eval application.
    \param arg1 is the input xml file, defaults to ../xml/test-dict.xml
*/
int main(int argc, char* argv[]) {
  char* xmlInput = "../xml/test-dict.xml";
  if (argc >= 2) {
    xmlInput = argv[1];
    if (!(strcmp(xmlInput, "-help"))) {  // instructions
      std::cout << 
      "First argument is xml file (id dictionary) to be parsed; defaults to test-dict.xml" 
                << std::endl;
      exit(0);
    }
  }

  xml::XmlParser* parser = new xml::XmlParser();
  DOM_Document doc = parser->parse(xmlInput);
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
  std::cout << "Dictionary " << dict->getDictName() << " made" << std::endl;

  
  // Attempt to verify it
  bool valid = dict->isValid();
  std::cout << "Dictionary is " << ( (valid) ? "VALID" : "INVALID" );
  std::cout << std::endl;

  // Test query functions
  testQuery(dict);

  // Register dictionary
  xmlUtil::IdDictMan *man = xmlUtil::IdDictMan::getPointer();
  xmlUtil::IdDictMan::RetCode ret = man->registerDict(dict);

  switch(ret) {
  case xmlUtil::IdDictMan::SUCCESS: 
    std::cout << "Dictionary successfully registered" << std::endl;
    break;
  case xmlUtil::IdDictMan::DUPLICATE: 
    std::cout << "Duplicate dictionary registration" << std::endl;
    break;
  default: 
    std::cout << "Dictionary registration failed" << std::endl;
    exit(ret);
  }

  // Check that we can find it again
  xmlUtil::IdDict *dictAgain = man->findDict(dict->getDictName());
  std::cout << ((dictAgain == dict) ? "Found dict" : "Failed to find dict") <<
    std::endl;

  // Look for IdConverter
  DOM_Element converterElt = 
    xml::Dom::findFirstChildByName(docElt, "idConverter");
  if (converterElt == DOM_Element()) {
    std::cout << "No id converter found " << std::endl;
    exit(0);
  }
  std::cout << "Id converter found" << std::endl;

  // Build representation
  xmlUtil::IdConverter *converter = new xmlUtil::IdConverter(converterElt);

  // What's it got?
  converter->displayConversions(std::cout);

  // Check it out
  return (testConverter(converter));
}

void testQuery(xmlUtil::IdDict* dict) {

  using namespace xmlUtil;  // for diagnostic output operators

  xmlUtil::Identifier okId(3);
  xmlUtil::Identifier badId(4);

  okId[0] = 0;
  okId[1] = 1;
  okId[2] = 1;

  badId[0] = 0;
  badId[1] = 0;
  badId[2] = 27;
  badId[3] = 1;

  bool valid = dict->idOk(okId);
  std::cout << "okId is " <<  ( (valid) ? "VALID" : "INVALID" );
  std::cout << std::endl;

  std::cout << "Original: " << okId << std::endl;
  std::cout << "Reconstituted: " << std::endl << (*(dict->getNamedId(okId)));

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

  std::cout << "Original named Id: " << nIdOk;
  Identifier *stripped = nIdOk.stripNames();
  std::cout << "..stripped: " << (*stripped);
  std::cout << "..reconstitued: " << std::endl 
            << (*(dict->getNamedId(*stripped)));

  valid = dict->idOk(nIdBadName);
  std::cout << "nIdBadName is " <<  ( (valid) ? "VALID" : "INVALID" );
  std::cout << std::endl;
}

int testConverter(xmlUtil::IdConverter* converter) {

  using namespace xmlUtil;  // to get diagnostic output functions

  if (!(converter->isConsistent()) ) {
    std::cout << "Converter is inconsistent" << std::endl;
    return (0);
  }
  std::cout << "Converter OK" << std::endl;

  // Check out some particular Identifier and NamedId objects
  xmlUtil::NamedId nId;
  nId.addField("fTopObjects", 0);   // LAT
  nId.addField("fLATObjects", 0);   // Towers
  nId.addField("fTowerY", 3);

  std::cout << "Before conversion:" << std::endl << nId;

  std::cout << "After conversion (shouldn't change):" << std::endl;
  
  std::cout << (*(converter->convert(&nId)));

  nId.addField("fTowerX", 2);
  nId.addField("fTowerObjects", 0);    // CAL
  nId.addField("fLayer", 3);
  nId.addField("fMeasure", 0);
  nId.addField("fCALLog", 4);
  nId.addField("fCALLogCmp", 2);

  std::cout << "Before conversion:" << std::endl << nId;

  std::cout << "After conversion (should truncate):" << std::endl;
  
  std::cout << (*(converter->convert(&nId)));

  nId.popField(5);
  nId.addField("fTowerObjects", 1);   //TKR
  nId.addField("fTKRTray", 6);
  nId.addField("fMeasure", 1);
  nId.addField("fTKRTrayCmp", 3);

  std::cout << "Before conversion:" << std::endl << nId;

  std::cout << "After conversion (should compress):" << std::endl;
  
  std::cout << (*(converter->convert(&nId)));

  nId.popField(7);
  nId.addField("fLATObjects", 1);   // ACD
  nId.addField("fACDFace", 3);
  nId.addField("fRow", 2);
  nId.addField("fCol", 3);

  std::cout << "Before conversion:" << std::endl << nId;

  std::cout << "After conversion (should disappear):" << std::endl;
  
  std::cout << (*(converter->convert(&nId)));

  // Try last one again, but on Identifier rather than named id
  xmlUtil::Identifier *id = nId.stripNames();

  std::cout << "Before conversion (Identifier):" << std::endl << (*id);
  std::cout << "After conversion (should disappear):" << std::endl;
  std::cout << (*(converter->convert(id)));

  return 1;
}


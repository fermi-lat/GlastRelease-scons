// $Header $
/*! \file See if DocMan and GDDDocMan perform properly.  Make a
    trivial DocClient which just outputs an element's tagname
    and the value of the  "substituted" attribute, if there is
    one.  Then parse document twice, once using DocMan and once
    using GDDDocMan.
*/
#include "xmlUtil/docMan/GDDDocMan.h"
#include "xml/docMan/DocClient.h"
#include <xercesc/dom/DOMElement.hpp>
//#include <xercesc/dom/DOMString.hpp>
#include "xml/Dom.h"

#include <string>
#include <iostream>
using namespace xmlUtil;
XERCES_CPP_NAMESPACE_USE

class MiniClient : public xml::DocClient {
public:
  MiniClient() {m_name = "miniClient";}
  MiniClient(std::string name) : m_name(name) {}
  const std::string& getName() { return m_name;}
  void handleChild(DOMNode* node);
  ~MiniClient() {};
private:
  std::string m_name;
};

void MiniClient::handleChild(DOMNode* node) {
  if (node->getNodeType() == DOMNode::ELEMENT_NODE) {
    DOMElement* elt = static_cast<DOMElement*>(node);
    //    DOMString eltName = node.getNodeName();
    std::string eltName = xml::Dom::getTagName(elt);
    std::cout << m_name 
              << " found element <" << eltName << ">" 
              << std::endl;
    std::string att = xml::Dom::getAttribute(elt,"substituted");
    if (att.size() > 0) {
      std::cout << "substituted = " << att << std::endl;
    }
  }
  else         {// tilt!  
    std::cout << "User client passed non-element node:  TILT " << std::endl;
  }
}

int main(int argc, char* argv[]) {


  std::string defFile("../xml/flight.xml");
  std::string* pFileName = &defFile;
  if (argc >= 2) pFileName = new std::string(argv[1]);

  // Need a client
  MiniClient* myClient = new MiniClient("firstClient");
  MiniClient* myClient2 = new MiniClient("secondClient");

  // First the standard DocMan
  std::cout << "Output from DocMan parsing " << std::endl;
  xml::DocMan* pMan = xml::DocMan::getPointer();
  pMan->regClient("constants", myClient);
  pMan->regClient("section", myClient);
  pMan->regClient("section", myClient2);
  pMan->regClient("source", myClient);  
  pMan->parse(*pFileName, "gdd");

  // Now GDDDocMan
  std::cout << std::endl << "Output from GGDDocMan parsing " << std::endl;
  xmlUtil::GDDDocMan* pGDDMan = xmlUtil::GDDDocMan::getPointer();
  pGDDMan->regClient("constants", myClient);
  pGDDMan->regClient("section", myClient);
  pGDDMan->regClient("idDict", myClient);
  pGDDMan->parse(*pFileName);
    
  return 0;
}

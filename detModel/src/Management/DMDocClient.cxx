#include "detModel/Management/DMDocClient.h"
#include "detModel/Utilities/Global.h"

#include <xercesc/dom/DOMElement.hpp>
#include "xmlBase/Dom.h"

#include <string>
#include <iostream>

namespace detModel{
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMNode;
  void DMDocClient::handleChild(DOMNode* node) {
    
    std::string name = xmlBase::Dom::getNodeName(node);

    if(name == "constants")
      {
	m_constants = node;
	std::cout << "detModel Doc Client: registered the constants!" << std::endl; 
      }
    else if (name == "section")
      {
	m_sections = node;
	std::cout << "detModel Doc Client: registered the sections!" << std::endl; 
      }
    else if (name == "materials")
      {
	m_materials = node;
	std::cout << "detModel Doc Client: registered the materials!" << std::endl; 
      }
    else if (name == "idDict")
      {
	m_idDict = node;
	std::cout << 
       "detModel Doc Client: registered the id dictionary!" << std::endl; 
      }
    else
      {
	detAbort ("detModel Doc Client asked for an unknown DOM element!");
      }
  }
}



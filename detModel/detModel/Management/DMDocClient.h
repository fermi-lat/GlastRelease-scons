#ifndef DMDOCCLIENT_H
#define DMDOCCLIENT_H
#include <string>
#include <iostream>

#include "xml/docMan/DocClient.h"
#include <xercesc/dom/DOMNode.hpp>
// #include <xercesc/dom/DOMString.hpp>
#include "xml/Dom.h"

namespace detModel{

  /*
    This is a concrete DocClient to be used by detModel
    to access xml files for parsing
  */

  class DMDocClient : public xml::DocClient {
  public:
    DMDocClient(){m_name = "detModelClient";}
    DMDocClient(std::string name) : m_name(name) {}
    
    const std::string& getName() { return m_name;}
    
    void handleChild(XERCES_CPP_NAMESPACE_QUALIFIER DOMNode* node);
    
    //    const DOM_Node* getConstants() {return &m_constants;}
    //    const DOM_Node* getSections() {return &m_sections;}
    //    const DOM_Node* getMaterials() {return &m_materials;}
    const XERCES_CPP_NAMESPACE_QUALIFIER DOMNode* getConstants() 
    {return m_constants;}
    const XERCES_CPP_NAMESPACE_QUALIFIER DOMNode* getSections() 
    {return m_sections;}
    const XERCES_CPP_NAMESPACE_QUALIFIER DOMNode* getMaterials() 
    {return m_materials;}
    XERCES_CPP_NAMESPACE_QUALIFIER DOMNode* getIdDictionary() 
    {return m_idDict;}

    ~DMDocClient() {};

  private:
    /// This is the name of the client
    std::string m_name;
    /// This is the DOM_Node for constants
    XERCES_CPP_NAMESPACE_QUALIFIER DOMNode* m_constants;
    /// This is the DOM_Node for sections
    XERCES_CPP_NAMESPACE_QUALIFIER DOMNode* m_sections;
    /// This is the DOM_Node for materials
    XERCES_CPP_NAMESPACE_QUALIFIER DOMNode* m_materials;
    /// This is the DOM_Node for id dictionary
    XERCES_CPP_NAMESPACE_QUALIFIER DOMNode* m_idDict;
  };
  
}
#endif

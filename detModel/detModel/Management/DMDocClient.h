#ifndef DMDOCCLIENT_H
#define DMDOCCLIENT_H
#include <string>
#include <iostream>

#include "xml/docMan/DocClient.h"
#include <xercesc/dom/DOM_Node.hpp>
#include <xercesc/dom/DOMString.hpp>
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
    
    void handleChild(DOM_Node node);
    
    //    const DOM_Node* getConstants() {return &m_constants;}
    //    const DOM_Node* getSections() {return &m_sections;}
    //    const DOM_Node* getMaterials() {return &m_materials;}
    const DOM_Node getConstants() {return m_constants;}
    const DOM_Node getSections() {return m_sections;}
    const DOM_Node getMaterials() {return m_materials;}
    const DOM_Node getIdDictionary() {return m_idDict;}

    ~DMDocClient() {};

  private:
    /// This is the name of the client
    std::string m_name;
    /// This is the DOM_Node for constants
    DOM_Node m_constants;
    /// This is the DOM_Node for sections
    DOM_Node m_sections;
    /// This is the DOM_Node for materials
    DOM_Node m_materials;
    /// This is the DOM_Node for id dictionary
    DOM_Node m_idDict;
  };
  
}
#endif

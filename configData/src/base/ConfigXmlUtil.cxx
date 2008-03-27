#include "configData/base/ConfigXmlUtil.h"

#include "xmlBase/Dom.h"
#include "xmlBase/XmlParser.h"

#include "xercesc/dom/DOMElement.hpp"
#include "xercesc/dom/DOMDocument.hpp"
#include "xercesc/dom/DOMImplementation.hpp"
#include "xercesc/util/XercesDefs.hpp"
#include <xercesc/util/XMLString.hpp>


namespace configData {

  namespace XmlUtil {
 
    void init() {
      static xmlBase::XmlParser* parser(0);
      if ( parser == 0 ) {
	parser = new xmlBase::XmlParser;
      }
    }

    DOMElement* makeDocument(const char* name) {

      using XERCES_CPP_NAMESPACE_QUALIFIER DOMImplementation;
      using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
      using XERCES_CPP_NAMESPACE_QUALIFIER XMLString;
      
      XMLCh tempStr[100];
      XMLString::transcode(name,tempStr,99);
      DOMElement* newChild = DOMImplementation::getImplementation()->createDocument()->createElement(tempStr);
      newChild->getOwnerDocument()->appendChild(newChild);
      return newChild;
    }


    DOMElement* makeChildNode(DOMElement& domNode, const char* name) {
      
      using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
      using XERCES_CPP_NAMESPACE_QUALIFIER XMLString;
      
      XMLCh tempStr[100];
      XMLString::transcode(name,tempStr,99);
      DOMElement* newChild = domNode.getOwnerDocument()->createElement(tempStr);
      domNode.appendChild(newChild);
      return newChild;
    }

    DOMElement* makeChildNodeWithContent(DOMElement& domNode, const char* name, const char* content) {
      DOMElement* cNode = makeChildNode(domNode,name);
      if ( cNode == 0 ) return 0;
      using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
      using XERCES_CPP_NAMESPACE_QUALIFIER XMLString;      
      XMLCh tempStr[100];
      XMLString::transcode(content,tempStr,99);
      cNode->setTextContent(tempStr);
      return cNode;
    }

    DOMElement* makeChildNodeWithContent(DOMElement& domNode, const char* name, unsigned int content) {
      char tmp[20];
      sprintf(tmp,"%d",content);	
      return makeChildNodeWithContent(domNode,name,tmp);
    }
    
    DOMElement* makeChildNodeWithContent(DOMElement& domNode, const char* name, void* content) {
      char tmp[20];
      sprintf(tmp,"%0x",(unsigned int)content);	
      return makeChildNodeWithContent(domNode,name,tmp);
    }

    DOMElement* findChildByName(DOMElement& parent, const char* name ) {
      return xmlBase::Dom::findFirstChildByName(&parent,name);
    }
    
    bool getAttribute(DOMElement& domNode, const char* attrName, std::string& theVal) {
      try { theVal = xmlBase::Dom::getAttribute(&domNode,attrName); }
      catch ( std::exception& ) { return false; }
      return true;
    }

    void getText(DOMElement& domNode, std::string& theVal) {
      theVal = xmlBase::Dom::getText(&domNode);
    }

    bool writeIt(DOMElement& doc,const char* fileName) {
      return xmlBase::Dom::writeIt(doc.getOwnerDocument(),fileName);
    }

  }

}



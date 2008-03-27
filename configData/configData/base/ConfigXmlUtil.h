#ifndef ConfigXmlUtil_h
#define ConfigXmlUtil_h 

// Include xerces
#include <xercesc/util/XercesDefs.hpp>

XERCES_CPP_NAMESPACE_BEGIN
class  DOMElement;
XERCES_CPP_NAMESPACE_END

using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;

/** 
 *
 * @brief A few functions to help manage xml
 *
 */

namespace configData {

  namespace XmlUtil {
    
    void init();
    
    DOMElement* makeDocument(const char* topNode);

    DOMElement* makeChildNode(DOMElement& domNode, const char* name);

    DOMElement* makeChildNodeWithContent(DOMElement& domNode, const char* name, const char* content);

    DOMElement* makeChildNodeWithContent(DOMElement& domNode, const char* name, unsigned int content);
    
    DOMElement* makeChildNodeWithContent(DOMElement& domNode, const char* name, void* content);
   
    DOMElement* findChildByName(DOMElement& parent, const char* name );
    
   
    void addAttribute(DOMElement& domNode, const char* name, const char* val);
    void addAttribute(DOMElement& domNode, const char* name, const int& val);
    void addAttribute(DOMElement& domNode, const char* name, const double& val);
    
    bool writeIt(DOMElement& doc, const char* fileName);
  }

};

#endif

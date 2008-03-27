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

    // Make an XmlParser, to initialize xerces
    void init();
    
    // Make a DOMDocument
    DOMElement* makeDocument(const char* topNode);

    // Make a DOMElement and attach it as a child to another element
    DOMElement* makeChildNode(DOMElement& domNode, const char* name);

    // Make a DOMElement and attach it as a child to another element
    // Also, set the text content
    DOMElement* makeChildNodeWithContent(DOMElement& domNode, const char* name, const char* content);

    // Make a DOMElement and attach it as a child to another element
    // Also, set the text content    
    DOMElement* makeChildNodeWithContent(DOMElement& domNode, const char* name, unsigned int content);

    // Make a DOMElement and attach it as a child to another element
    // Also, set the text content        
    DOMElement* makeChildNodeWithContent(DOMElement& domNode, const char* name, void* content);
   
    // Find a DOMElement with a given name in the set of children of a given element
    DOMElement* findChildByName(DOMElement& parent, const char* name );
    
    // Write a DOMDocument to a file
    bool writeIt(DOMElement& doc, const char* fileName);
  }

};

#endif

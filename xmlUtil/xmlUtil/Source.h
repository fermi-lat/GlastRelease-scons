// $Header$
// Author: J. Bogart

#ifndef XMLUTIL_SOURCE_H
#define XMLUTIL_SOURCE_H

#include <string>
// #include <xercesc/dom/DOM_Element.hpp>
// #include <xercesc/dom/DOM_Document.hpp>

#include <xercesc/util/XercesDefs.hpp>

XERCES_CPP_NAMESPACE_BEGIN
class DOMElement;
class DOMDocument;
XERCES_CPP_NAMESPACE_END

namespace xmlUtil {
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;

  /*! This class provides an interface to the <source> element of
       gdd.dtd.  For now only provide a way to write one.  Someday
       might also provide interface to read pre-existing ones.
  */
  class Source {
  public:
    Source(DOMDocument* doc, const char * creator, 
           const char* creatorCVSID);
    ~Source() {};

     void add(DOMElement* after = 0);

  private:
    DOMDocument* m_doc;
    DOMElement*  m_me;  // the element as represented in DOM
    std::string  m_creator;
    std::string  m_myCVSID;
  };
}
#endif

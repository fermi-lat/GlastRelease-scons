// $Header$
// Author: J. Bogart

#ifndef XMLUTIL_SOURCE_H
#define XMLUTIL_SOURCE_H

#include <string>
#include <dom/DOM_Element.hpp>
#include <dom/DOM_Document.hpp>

namespace xmlUtil {
  /*! This class provides an interface to the <source> element of
       gdd.dtd.  For now only provide a way to write one.  Someday
       might also provide interface to read pre-existing ones.
  */
  class Source {
  public:
    Source(DOM_Document doc, const char * creator, 
           const char* creatorCVSID);
    ~Source() {};

     void add(const DOM_Element after = DOM_Element());

  private:
    DOM_Document m_doc;
    DOM_Element  m_me;  // the element as represented in DOM
    std::string  m_creator;
    std::string  m_myCVSID;
  };
}
#endif

// $Header$
#ifndef XMLUTIL_DOCCLIENT_H
#define XMLUTIL_DOCCLIENT_H

#include "dom/DOM_Element.hpp"

namespace xmlUtil {
  //! Pure virtual class whose only purpose is to provide a 
  //! handler for clients of DocMan
  class DocClient {
  public:
    DocClient() {}
    virtual ~DocClient() {}

    //! Must be overridden by derived classes
    virtual void handleChild(DOM_Node node)=0;  

    //! Give DocMan a way to identify clients and detect duplicates; 
    //! must be overridden by derived classes
    virtual const std::string& getName()=0;
  };
}

#endif

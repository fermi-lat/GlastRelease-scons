// $Header$
#ifndef XMLUTIL_GDDDOCMAN_H
#define XMLUTIL_GDDDOCMAN_H

#include "xmlUtil/docMan/DocMan.h"
#include "xmlUtil/docMan/DocClient.h"
#include "xmlUtil/Substitute.h"

namespace xmlUtil {
  //! GDDDocMan knows about constants sections.  As part of construction
  //! it signs up an object it owns as a handler for <constants> so it 
  //! can do evaluation and substitution before any other clients see 
  //! the DOM
  class GDDDocMan : public DocMan {
  public:
    static  GDDDocMan* getPointer();

    // Allow caller option of saving (or not) notes associated with constants 
    virtual bool parse(const std::string& filename, bool saveNotes = true);

    virtual ~GDDDocMan() {delete m_constsClient;}

  protected:
    GDDDocMan::GDDDocMan();            // TO BE WRITTEN
    
    //! Nested class to do constansts handling; inherits from DocClient
    class ConstsClient : public DocClient {
    public:
      virtual const std::string& getName() {
        return name;
      }
      ConstsClient(GDDDocMan *owner) : m_owner(owner), m_sub(0) {
        name=std::string("ConstsEvalSub");
        m_doc = 0;}
      virtual ~ConstsClient() {}

      void handleChild(DOM_Node elt);       /* TO BE WRITTEN */
      bool m_saveNotes;
    private:
      std::string name;
      GDDDocMan* m_owner;
      DOM_Document m_doc;
      Substitute* m_sub;
    };
    static GDDDocMan* m_me;
    ConstsClient* m_constsClient;
  };
}
#endif

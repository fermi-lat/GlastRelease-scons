
// $Header$
#ifndef XMLUTIL_GDDDOCMAN_H
#define XMLUTIL_GDDDOCMAN_H

#include "xmlBase/docMan/DocMan.h"
#include "xmlBase/docMan/DocClient.h"
#include "xmlUtil/Substitute.h"

namespace xmlUtil {
  //! GDDDocMan knows about constants sections.  As part of construction
  //! it signs up an object it owns as a handler for <constants> so it 
  //! can do evaluation and substitution before any other clients see 
  //! the DOM
  class GDDDocMan : public xmlBase::DocMan {
  public:
    static  GDDDocMan* getPointer();

    // Allow caller option of saving (or not) notes associated with constants 
    // docType should be gdd (should we even allow anything else?)
    virtual bool parse(const std::string& filename, 
                       bool saveNotes = true,
                       const std::string& docType=std::string("gdd"));
    const std::string& getCVSid() const;
    const std::string& getDTDversion() const;

    virtual ~GDDDocMan() {delete m_constsClient;}

  protected:
    GDDDocMan::GDDDocMan();            // TO BE WRITTEN
    
    //! Nested class to do constansts handling; inherits from DocClient
    class ConstsClient : public xmlBase::DocClient {
    public:
      virtual const std::string& getName() {
        return name;
      }
      ConstsClient(GDDDocMan *owner) : m_owner(owner), m_sub(0),
                                       CVSid(std::string("")),
                                       DTDversion(std::string(""))
      {
        name=std::string("ConstsEvalSub");
        m_doc = 0;
      }

      virtual ~ConstsClient() {}

      void handleChild(XERCES_CPP_NAMESPACE_QUALIFIER DOMNode* node); 
      bool m_saveNotes;
      const std::string& getCVSid() const {return CVSid;}
      const std::string& getDTDversion() const {return DTDversion;}
    private:
      std::string name;
      GDDDocMan* m_owner;
      DOMDocument* m_doc;
      Substitute* m_sub;
      std::string CVSid;
      std::string DTDversion;
    };
    static GDDDocMan* m_me;
    ConstsClient* m_constsClient;
  };
}
#endif

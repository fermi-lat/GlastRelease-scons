// $Header$
#ifndef XMLUTIL_DICTVALIDVISITOR_H
#define XMLUTIL_DICTVALIDVISITOR_H

#include "xmlUtil/id/DictVisitor.h"

namespace xmlUtil {

  class DictValidVisitor : public DictVisitor {
  public:
    DictValidVisitor() : m_dict(0), m_state(eReady) {}
    bool         visitDict();

    virtual bool visitDict(IdDict *dict);
    virtual bool visitNode(DictNode *node);

    virtual bool visitField(DictField *) {return true;}
    virtual bool visitFieldMan(DictFieldMan *) {return true;}

    bool wasValid() {return (m_state == eValid); }
    
  private:
    enum State {eReady, eValid, eInvalid};
    IdDict* m_dict;
    State   m_state;
    
  };
}

#endif

#ifndef TrgConfigParser_h
#define TrgConfigParser_h

#include <xercesc/dom/DOMElement.hpp>
#include "xmlBase/Dom.h"

#include "TrgConfig.h"
#include <string>
#include <list>
#include <map>
using namespace xmlBase;

class TrgConfigParser{
 public:
  TrgConfigParser(const char* filename);
  TrgConfigParser(DOMElement* topElt);
  int parse(TrgConfig* tcf);
 private:
  const std::string removeWhitespace(const std::string) const;
  unsigned long content(DOMElement* el,const char* tag);
  DOMElement* m_topElt;
};

#endif
   

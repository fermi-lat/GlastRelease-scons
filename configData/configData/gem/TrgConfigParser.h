/** @class TrgConfigParser
    @author M. Kocian

    Class to parse a latc xml file for GEM information to fill a TrgConfig class object.
*/

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
  /// Default constructor
  TrgConfigParser();
  /// Constructor which takes the latc file name as a parameter
  TrgConfigParser(const char* filename);
  /// Constructor which takes the top DOMElement as a parameter
  TrgConfigParser(DOMElement* topElt,bool allowMissing=false);
  /// Function to fill TrgConfig object *tcf
  int parse(TrgConfig* tcf);
  /// Function to fill TrgConfig with contents of xml file filename
  int parse(TrgConfig* tcf, const char* filename);
  /// Function to fill TrgConfig with contents of DOMElement el;
  int parse(TrgConfig* tcf, DOMElement* el);
 private:
  /// Remove spaces from xml strings
  const std::string removeWhitespace(const std::string) const;
  /// return the text content of tag tag in DOMElement el
  unsigned long content(DOMElement* el,const char* tag);
  /// split the tag for ROI into individual tile tags
  const std::string maketag(const char*,const int)const;
  DOMElement* m_topElt;
  bool m_allowMissing;
};

#endif
   

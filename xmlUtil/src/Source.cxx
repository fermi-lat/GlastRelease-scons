// $Header$

#include "xmlUtil/Source.h"
// #include <xercesc/dom/DOMString.hpp>
#include <xercesc/dom/DOMText.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/util/XMLString.hpp>
#include "xmlBase/Dom.h"


namespace { 
  // private utility to strip off $ 
  char * stripDollar(char* toStrip) {
    const char dollar = '$';
    
    char   nothing = 0;
    unsigned len = strlen(toStrip);
    
    if (dollar == toStrip[len - 1]) {
      len--;
      toStrip[len] = nothing;
    }
    if (dollar == toStrip[0]) {
      return (toStrip + 1);
    }
    return toStrip;
  }
}      // end of blank namespace

namespace xmlUtil {
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMText;
  using XERCES_CPP_NAMESPACE_QUALIFIER XMLString;
  Source::Source(DOMDocument* doc, const char * creator,
      const char* creatorCVSID) : m_doc(doc), m_creator(creator)  {
    
    char * tmp = new char[strlen(creatorCVSID) + 1];
    strcpy(tmp, creatorCVSID);
    m_myCVSID = std::string(stripDollar(tmp));
    delete [] tmp;
  }

  void Source::add(DOMElement* parentArg) {
    static char * unknownId = "unidentified input file";
    DOMElement* docElt = m_doc->getDocumentElement();
    DOMElement* parent = parentArg;
    if (parent == 0) {
      parent = docElt;
    }

    XMLCh* xmlchSource = XMLString::transcode("source");
    DOMElement* source = m_doc->createElement(xmlchSource);
    XMLString::release(&xmlchSource);

    // Document element should have a CVSid attribute
    std::string idAtt = xmlBase::Dom::getAttribute(docElt, "CVSid");

    char * strippedRaw;
    char * myRaw = 0;

    unsigned int len = idAtt.size();
    if (len > 0 ) {
      const char * raw = idAtt.c_str();
      myRaw = new char[len + 1];
      strcpy(myRaw, raw);
      strippedRaw = stripDollar(myRaw);
    }
    else strippedRaw = unknownId;

    xmlBase::Dom::addAttribute(source, "inputId", std::string(strippedRaw));

    // Already have creator ID as std::string
    xmlBase::Dom::addAttribute(source, "creatorId", m_myCVSID);
    
    std::string theText = 
      std::string("Do not edit! This file automatically created by ");
    theText.append(m_creator);
    XMLCh* xmlchText = XMLString::transcode(theText.c_str());
    
    DOMText*    sourceText = m_doc->createTextNode(xmlchText);
    XMLString::release(&xmlchText);

    source->appendChild(sourceText);

    if (myRaw) delete [] myRaw;
    // Finally need to insert this node as the first child of specified elt
    parent->insertBefore(source, parent->getFirstChild());
  }
}

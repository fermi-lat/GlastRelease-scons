// $Header$
// Author: J. Bogart

#ifndef XMLUTIL_CONSTANTS_H
#define XMLUTIL_CONSTANTS_H

/*
#include <xercesc/dom/DOM_Element.hpp>
#include <xercesc/dom/DOM_Document.hpp>
*/
#include <xercesc/util/XercesDefs.hpp>
#include <vector>
XERCES_CPP_NAMESPACE_BEGIN
class DOMElement;
class DOMDocument;
class DOMTreeWalker;
XERCES_CPP_NAMESPACE_END

namespace xmlUtil {
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMElement;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument;
  using XERCES_CPP_NAMESPACE_QUALIFIER DOMTreeWalker;
  //! This class is a place to keep functions which manipulate the
  //! <constants> element and its children in the DOM representation.
  class Constants {
  public:
    Constants(DOMDocument* doc);
    ~Constants();
    //! Fix up and length units that need fixing in <prim> elements
    void normalizePrimary();
    //! After constants have been evaluated, remove child elements
    //! (except <notes> when argument says to keep them)
    void pruneConstants(const bool keepNotes = false);

    //! For each derCategory, evaluate each const child and save
    //! in "value" attribute
    void evalConstants();
  private:
    void handleEnergies(std::vector<DOMElement*> saved);

    DOMDocument* m_doc;          // document element
    DOMElement*  m_constants;    // reference to <constants> element
    DOMTreeWalker* m_walker;
  };
}
#endif

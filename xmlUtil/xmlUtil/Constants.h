// $Header$
// Author: J. Bogart

#ifndef XMLUTIL_CONSTANTS_H
#define XMLUTIL_CONSTANTS_H

#include <dom/DOM_Element.hpp>
#include <dom/DOM_Document.hpp>

namespace xmlUtil {
  //! This class is a place to keep functions which manipulate the
  //! <constants> element and its children in the DOM representation.
  class Constants {
  public:
    Constants(DOM_Document doc);
    //! Fix up and length units that need fixing in <prim> elements
    void normalizePrimary();
    //! After constants have been evaluated, remove child elements
    //! (except <notes> when argument says to keep them)
    void pruneConstants(const bool keepNotes = false);

    //! For each derCategory, evaluate each const child and save
    //! in "value" attribute
    void evalConstants();
  private:
    DOM_Document m_doc;          // document element
    DOM_Element  m_constants;    // reference to <constants> element
  };
}
#endif

// $Header$
#ifndef XMLUTIL_SUBSTITUTE_H 
#define XMLUTIL_SUBSTITUTE_H

#include <dom/DOM_Document.hpp>
#include <dom/DOMString.hpp>
#include <dom/DOM_TreeWalker.hpp>
#include <string>


namespace xmlUtil {
  /*! The Substitute class will look up references to values (in
      attributes of type IDREF in a given subtree of a document)
      and replace the attribute with another that has as value
      the value of the reference.  e.g.,
             ... XREF="CsiLength"   
      could get replaced by
             ... X = "33.12"
  */
  class Substitute {
  public: 

    /*!  A Substitute object is configured by specifying the document
         in which substitutions are to occur and the suffix which
         will be used to find attributes for which substitutions
         are requested.
         \param doc  The document for which substitutions are to be
         performed
         \param suffix The suffix which identitifies references in
         the source docoument (defaults to "REF")
    */
    Substitute::Substitute(DOM_Document doc, std::string suffix =
                           std::string("REF")); 

    /*! Do substitution for the subtree indicated, returning count.
        \b execute may be invoked repeatedly for different subtrees
        of the source document.
        \param treeTop  the root of the subtree for which substitutions
        are to be performed.
    */
    int execute(DOM_Element treeTop);

    //! Destructor.
    ~Substitute() {};

  private:
    DOM_Document m_doc;    /*!< Document on which substitutions may be done */
    int          m_count;  /*!< Keep track of \#subs done for one \b execute */
    std::string  m_suffix; /*!< Suffix identifying references */
    int          m_suffixLen; 
    int          m_notFound; /*!< count of id's not found-- shouldn't happen */

    /*! Handle a single element (invoked from public execute method) */
    void  sub(DOM_Element elt);
  };

}
#endif

// $Header$
#ifndef XMLUTIL_ARITH_H 
#define XMLUTIL_ARITH_H

#include <dom/DOM_Element.hpp>
#include <dom/DOMString.hpp>

namespace xmlUtil {
  //! The Arith class handles the xml arithmetic elements, allowing
  //! derived constants to be computed and their values put into
  //! the DOM representation of the XML file.
  class Arith {
  public: 
    //! The different element types handled (evaluated) by this class
    enum  eTagType {
      ETAG_const = 0,
      ETAG_refer, 
      ETAG_add,
      ETAG_minus,
      ETAG_mul,
      ETAG_quo,
      ETAG_uminus,
      ETAG_max,
      ETAG_n 
    } ;

    Arith::Arith(const DOM_Element elt);
    //! Return a double, the result of performing all the operations
    //! indicated by the element and its children. 
    double evaluate();

    //! add "value" attribute into DOM. 
    /*! This is a no-op for all but const elements. */
    void   saveValue(); 

    //! Return scale factor due to units specification.  
    /*! For now only applies to lengths, which should be converted 
        to mm (always returns 1.0 for other kinds of constants)
    */
    static double getUnitsScale(const DOM_Element& elt);
    ~Arith() {};

  private:
    static void init(); /*< Sets up static strings */
    double add();       /*< evaluate an \b add element */
    double minus();     /*< evaluate a \b minus element */
    double mul();       /*< evaluate a \b mul element */
    double quo();       /*< evalauate a \b quo element */
    double myMax();     /*< evaluate a \b max element */

    DOM_Element  m_elt; /*< element which may be evaluated */
    int     m_tag;      /*< element type (tagname) */
    double  m_number;   /*< value of element after evaluation */
    bool    m_evaluated;/*< flag to keep track of whether eval has been done*/

    typedef DOMString * ptrDOMString;/*< keep some static strings so we don't
                                        have to build them over and over */
    static  ptrDOMString *typeNames;
    static  ptrDOMString valString;   // a DOMString containing "value"
    static  ptrDOMString refToString; // ..and so forth
    static  ptrDOMString constString;
    static  ptrDOMString referString;
    static  ptrDOMString addString;
    static  ptrDOMString minusString;
    static  ptrDOMString uminusString;
    static  ptrDOMString mulString;
    static  ptrDOMString quoString;
    static  ptrDOMString maxString;
    static  ptrDOMString notesString;
    static  ptrDOMString lengthString;
    static  ptrDOMString cmString;
    static  ptrDOMString mString;
  };

}
#endif

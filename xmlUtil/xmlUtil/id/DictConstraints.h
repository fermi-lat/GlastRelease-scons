// $Header$

#ifndef XMLUTIL_DICTCONSTRAINTS_H
#define XMLUTIL_DICTCONSTRAINTS_H

#include <vector>
#include <dom/DOM_Element.hpp>
#include "xmlUtil/id/DictVisitor.h"
#include <set>

namespace xmlUtil {
  //! Represent constraints on values for an id dictionary node


  class DictConstraints {
  public: 

    //! Build constraints object from xml representation
    DictConstraints(DOM_Element elt);
    bool equals(const DictConstraints& other);
    ~DictConstraints(){if (m_valList) delete m_valList;}
  public:    
    //! Does the specified value satisfy the constraints?
    bool allowed(const unsigned value) const;

    //! Are the constraints of the argument a refinement (i.e.,
    //! define a subset of possible values)?
    bool allowed(const DictConstraints& other) const; 
    bool allowed(DictConstraints *other) const; 
    bool disjoint(const DictConstraints& other) const;      
    bool disjoint(DictConstraints* other) const;      

    unsigned getMin() const {return m_minVal;}
    unsigned getMax() const {return m_maxVal;}
    //    bool     isList() const {return (m_valList != 0);}

    //! Add our allowed values to specified set so we
    //! can use set operations to look for overlap
    void     insertValues(std::set<unsigned>& aSet) const;

    // DictNode needs to make copies of DictConstraints objects in 
    // some of its private constructors
    DictConstraints(const DictConstraints& toCopy);
    DictConstraints& operator=(const DictConstraints&);

  private:
    typedef std::vector<unsigned> DictValList;

    //! Utility used by copy constructor, copy assignment operator
    void deepCopy(const DictConstraints& toCopy);
    
    enum eStyle {
      ESTYLE_uninit = 0,
      ESTYLE_single = 1,
      ESTYLE_interval = 2,
      ESTYLE_list = 3
    };
    eStyle    m_style;    

    DictValList *m_valList;
    unsigned    m_minVal; //! Will be initialized correctly even for list
    unsigned    m_maxVal; //! Will be initialized correctly even for list

    //! Default constructor, not normally invoked
    DictConstraints() : m_style(ESTYLE_uninit), m_valList(0) {};

    eStyle   getType() const {return m_style;}

    DictConstraints(const unsigned soleValue);
    DictConstraints(const unsigned min, const unsigned max);
    DictConstraints(const DictValList * const  list);
    
  };       // end DictConstraints class
}          // end xmlUtil namespace

#endif

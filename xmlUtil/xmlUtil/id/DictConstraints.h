// $Header$

#ifndef XMLUTIL_DICTCONSTRAINTS_H
#define XMLUTIL_DICTCONSTRAINTS_H

#include <string>
#include <vector>
#include "dom/DOM_Element.hpp"

namespace xmlUtil {
  //! Represent constraints on values for an id dictionary node
  typedef std::vector<unsigned> DictValList;
  class DictConstraints {
  public: 
    //! Build constraints object from xml representation
    DictConstraints(DOM_Element elt);
    ~DictConstraints(){if (m_valList) delete m_valList;};

    //! Does the specified value satisfy the constraints?
    bool allowed(const unsigned value);

  private:
    // Not sure we need these
    friend class DictNode;
    friend class IdDict;

    // DictNode needs to make copies of DictConstraints objects in 
    // some of its private constructors
    DictConstraints(DictConstraints& toCopy);
    DictConstraints& operator=(const DictConstraints&);

    // Utility used by above
    void deepCopy(DictConstraints& toCopy);
    
    //! Default constructor, not normally invoked
    DictConstraints() : m_style(ESTYLE_uninit) {};

    unsigned getMin() {return m_minVal;}
    unsigned getMax() {return m_maxVal;}
    eStyle   getType() {return m_style;}

    DictConstraints(const unsigned soleValue);
    DictConstraints(const unsigned min, const unsigned max);
    DictConstraints(const DictValList const * list);
    unsigned         m_minVal, m_maxVal;
    DictValList *m_valList;
    enum eStyle {
      ESTYLE_uninit = 0,
      ESTYLE_single = 1,
      ESTYLE_interval = 2
      ESTYLE_list = 3;
    };
    eStyle    m_style;    
    
  };       // end DictConstraints class
}          // end xmlUtil namespace





#endif

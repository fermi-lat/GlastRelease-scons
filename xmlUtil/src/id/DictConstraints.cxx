// $Header$

#include <string>
#include <algorithm>
#include "xmlBase/Dom.h"
#include "xmlUtil/id/DictConstraints.h"
#include <xercesc/dom/DOMElement.hpp>


namespace xmlUtil {
  XERCES_CPP_NAMESPACE_USE
  /*! One and only public constructor.  Note that constraints specified
      as a list in xml, but describable in terms of a min and max
      bound will be represented in terms of the bounds.  The list
      representation is only used when there is no other option. */
  DictConstraints::DictConstraints(DOMElement* elt) : 
    m_style(ESTYLE_uninit), m_valList(0), m_minVal(0), m_maxVal(0)
  {
    using xmlBase::Dom;
    int   minVal, maxVal;

    std::string vType = Dom::getTagName(elt);
    if (vType == std::string("vEnumVal")) {
      m_minVal = minVal = Dom::getIntAttribute(elt, "value");
      m_maxVal = maxVal = minVal;
      m_style = ESTYLE_single;
    }
    else if (vType == std::string("vMinMax") ) {
      m_minVal = minVal = Dom::getIntAttribute(elt, "min");
      m_maxVal = maxVal = Dom::getIntAttribute(elt, "max");
      m_style = (m_minVal == m_maxVal) ? ESTYLE_single : ESTYLE_interval;
    }
    // Probably should revise this to sort the list before storing
    // For now, just update min and max correctly
    else if (vType == std::string("vList")) { 
      minVal = 0x7fffffff;
      maxVal = -1;

      std::vector<DOMElement*> children;
      // a vList should have only vEnumVal children (according to dtd)
      Dom::getChildrenByTagName(elt, "vEnumVal", children);
      unsigned nChild = children.size();
      if (nChild == 0) {
        m_style = ESTYLE_uninit;
        return;
      }
      m_valList = new std::vector<unsigned>(nChild);
      for (unsigned ix = 0; ix < nChild; ix++) {
        int temp;
        DOMElement*    childElt = children[ix];
        try {
          temp = Dom::getIntAttribute(childElt, "value");
          if (temp < 0) {
            delete m_valList;
            std::string msg =
              "From xmlUtil::DictConstraints::DictConstraints. vEnumVal has value < 0";
            throw xmlBase::WrongAttributeType(msg);
          }
          else {
            unsigned newVal = temp;
            m_valList->push_back(newVal);
            if (temp < minVal) minVal = temp;
            if (temp > maxVal) maxVal = temp;
          }
        }
        catch (xmlBase::DomException ex) {
          delete m_valList;
          std::cerr << "From xmlUtil::DictConstraints::DictConstraints " <<
            ex.getMsg();
          throw ex;
        }
      }     // end for loop
    }     // end val list case

    if ((minVal < 0 ) || (maxVal < minVal)) {
      if (m_valList) delete m_valList;
      std::string msg = "From xmlUtil::DictConstraints::DictConstraints ";
      msg += "illegal min and max ";
      throw xmlBase::WrongAttributeType(msg);
    }
    else {
      m_minVal = minVal;
      m_maxVal = maxVal;
    }
      
    // Check if it's really an interval after all.  Look for all
    // inbetween elements
    if (m_valList) {
      for (unsigned iy = m_minVal + 1; iy < m_maxVal; iy++) {
        // if there is anything we can't find, we really need a list
        if (std::find(m_valList->begin(), m_valList->end(), iy) == 
            m_valList->end()) {
          m_style = ESTYLE_list;
          return; 
        }
      }
      // Everything was there.  The bounds are sufficient.
      delete m_valList;
      m_style = (m_minVal == m_maxVal) ? ESTYLE_single : ESTYLE_interval;
      m_valList = 0;
    }
  }


  DictConstraints::DictConstraints(const unsigned soleValue) :
    m_style(ESTYLE_single), m_valList(0), m_minVal(soleValue), 
    m_maxVal(soleValue) {}

  DictConstraints::DictConstraints(const unsigned min, const unsigned max) :
    m_style(ESTYLE_interval), m_valList(0),
    m_minVal(min), m_maxVal(max) {}

  DictConstraints::DictConstraints(const DictValList * const list) :
    m_style(ESTYLE_list) {
    m_valList = new DictValList(*list);
    m_minVal = 0xffffffff;
    m_maxVal = 0;

    for (DictValList::const_iterator it = list->begin(); 
    it != list->end(); ++it) {
      const unsigned val = *it;
      if (val > m_maxVal) m_maxVal = val;
      if (val < m_minVal) m_minVal = val;
    }
  }

  DictConstraints::DictConstraints(const DictConstraints& toCopy) {
    deepCopy(toCopy);
  }

  void DictConstraints::deepCopy(const DictConstraints& toCopy) {
    // Nothing to it unless the DictValList is being used
    m_style = toCopy.m_style;
    m_valList = 0;
    m_minVal = toCopy.m_minVal;
    m_maxVal = toCopy.m_maxVal;
    if (toCopy.m_valList) {
      m_valList = new DictValList(*(toCopy.m_valList));
    }
  }


  DictConstraints& DictConstraints::operator=(const DictConstraints& d) {
    if (this != &d) {
      if (m_valList) delete m_valList;
      deepCopy(d);
    }
    return *this;
  }

  bool DictConstraints::equals(const DictConstraints& other) {
    if ((m_style != other.m_style) ||
        (m_minVal != other.m_minVal) ||
        (m_maxVal != other.m_maxVal) ) return false;

    if (m_style == ESTYLE_list) {
      return (*m_valList == *(other.m_valList));
    }
    else return true;
  }

  bool DictConstraints::allowed(const unsigned value) const {
    switch(m_style) {
    case ESTYLE_single: return (value == m_minVal);
    case ESTYLE_interval: 
      return ((value >= m_minVal) & (value <= m_maxVal));
    case ESTYLE_list:
      return 
        (std::find(m_valList->begin(), m_valList->end(), value) != 
         m_valList->end());
    default:
      
      return false;
    }
    return false;
  }
  
  bool DictConstraints::allowed(DictConstraints *other) const {
    return allowed(*other);
  }

  bool DictConstraints::allowed(const DictConstraints& other) const {

    bool boundsOk = ((other.m_minVal >= m_minVal) && 
                     (other.m_maxVal <= m_maxVal) );
    // if this isn't a list, this is all there is to check
    if (m_valList == 0) return boundsOk;
    // If either has a list of values, bounds must still be ok
    if (!(boundsOk)) return false;
                                                        
    // All that's left to check is case where this has a list
    // and  bounds are ok.  Just do the brute force thing and
    // check each possible value separately.  This would be horrible
    // if both were large enumerated sets, but it's very unlikely
    // that large sets will be described as a list rather than as an
    // interval

    else if (other.m_valList) {
      for (unsigned ix = 0; ix < (other.m_valList)->size(); ix++) {
        if (!allowed((*(other.m_valList))[ix])) return false;
      }
    } 
    else {  // other is interval-type constraint
      for (unsigned iy = other.m_minVal; iy <= other.m_maxVal; iy++) {
        if (!allowed(iy)) return false;
      }
    }
    return true;
  }

  bool DictConstraints::disjoint(const DictConstraints& other) const {
    if ((other.m_minVal > m_maxVal) ||
        (other.m_maxVal < m_minVal))  return true;
    // If either is an interval, we'll have an intersection
    if ((!m_valList) || (!other.m_valList)) return false;

    // Both are lists.  Check for duplicates one at a time
    for (DictValList::iterator it = m_valList->begin(); 
         it !=m_valList->end(); ++it) {
      unsigned val = *it;
      if (other.allowed(val)) return false;
    }
    return true;
  }

  bool DictConstraints::disjoint(DictConstraints *other) const {
    return disjoint(*other);
  }

  void DictConstraints::insertValues(std::set<unsigned>& aSet) const {
    if (m_valList == 0) {
      for (unsigned i = m_minVal; i <= m_maxVal; i++) {
        aSet.insert(i);
      }
    }
    else {
      aSet.insert(m_valList->begin(), m_valList->end());
    }
    return;
  }

}

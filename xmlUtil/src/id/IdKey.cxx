// $Header$

#include "xmlUtil/id/IdKey.h"
#include <cassert>

namespace {
  static const unsigned int maxVal = 254;
  static const unsigned int moveLeft = 256;
}

namespace xmlUtil {
  IdKey::IdKey(const std::vector<unsigned int> id) {
    fillFrom(id);
  }

  IdKey::IdKey(const std::deque<unsigned int> id) {
    fillFrom<std::deque<unsigned int> >(id);
  }

  // Compress each field to a byte.  Use the range 1-255 because
  // 0 is reserved; it means "absent".  If any field in the
  // input id has a value > 254 we're in trouble
  template<class Container>
  void IdKey::fillFrom(const Container& con) {

    typename Container::const_iterator start  = con.begin();
    typename Container::const_iterator stop  = con.end();

    const unsigned int inSize = con.size();

    m_key.reserve((inSize + 3) / 4);
    unsigned int out = 0;
    unsigned int iByte = 0;

    typename Container::const_iterator inIt = start;

    for (; inIt != stop; ++inIt) {
      unsigned int cur = *inIt + 1;

      if (cur > maxVal) {
        // not allowed
        assert(maxVal > cur);
      }
      out *= moveLeft;  
      out += cur; 
      iByte++; 
      if (iByte == 4) {  // store the word; re-init for next word
        m_key.push_back(out);
        out = 0;
        iByte = 0;
      }
    }

    if (iByte) {
      while (iByte < 4) {
        out *= moveLeft;
        iByte++;
      }
      m_key.push_back(out);
    }
  }

  // For sorting purposes, this must behave like < (less-than)
  // Internally keys are essentially vectors of unsigned ints.
  // Compare them component-wise.  Return 
  //     * when corresponding components differ       or
  //     * when one key is shorter than the other (shorter is less-than)
  bool IdKey::ltkey::operator()(IdKey key1, IdKey key2) const {

    IdKey::KeyIt it1 = key1.m_key.begin();
    IdKey::KeyIt it2 = key2.m_key.begin();

    while (it1 != key1.m_key.end()) {
      if (it2 == key2.m_key.end()) return false;
      if (*it1 != *it2) return (*it1 < *it2);
      ++it1; ++it2;
    }

    // If we got here the first key has been exhausted.
    return (it2 != key2.m_key.end());
  }

  IdKey::~IdKey() {m_key.clear(); m_key.resize(0);}
}
  
    

// $Header$

#include "CalibData/Cal/Xpos.h"

namespace CalibData {

  Xpos::Xpos(std::vector<float>* vals) {
    if (vals) {
      m_xpos.clear();
      m_xpos.reserve(vals->size());
      for (unsigned int i = 0; i < vals->size(); i++) {
        float val = (*vals)[i];
        m_xpos.push_back(val);
      }
    }
  }
  
  void Xpos::update(const Xpos* other) {
    m_xpos.clear();
    unsigned size = other->m_xpos.size();
    m_xpos.reserve(size);
    for (unsigned i = 0; i < size; i++) {
      m_xpos.push_back((other->m_xpos)[i]);
    }
  }

}

// $Header$

#include "CalibData/Cal/Xpos.h"

namespace CalibData {

  Xpos::Xpos(std::vector<float>* vals) {
    if (vals) {
      m_vals.clear();
      m_vals.reserve(vals->size());
      for (unsigned int i = 0; i < vals->size(); i++) {
        float val = (*vals)[i];
        m_vals.push_back(val);
      }
    }
  }
  
  void Xpos::update(const Xpos* other) {
    m_vals.clear();
    unsigned size = other->m_vals.size();
    m_vals.reserve(size);
    for (unsigned i = 0; i < size; i++) {
      m_vals.push_back((other->m_vals)[i]);
    }
  }

}

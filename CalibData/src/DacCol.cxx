// $Header$

#include "CalibData/DacCol.h"

namespace CalibData {

  DacCol::DacCol(std::vector<unsigned>* vals=0) {
    if (vals) {
      m_dacs.clear();
      for (unsigned int iDac = 0; iDac < vals->size(); iDac++) {
        m_dacs.push_back((*vals)[iDac]);
      }
    }
  }
  DacCol::DacCol(std::vector<int>* vals=0) {
    if (vals) {
      m_dacs.clear();
      for (unsigned int iDac = 0; iDac < vals->size(); iDac++) {
        unsigned val = (*vals)[iDac];
        m_dacs.push_back(val);
      }
    }
  }
  
  void DacCol::update(const DacCol* other) {
    m_dacs.clear();
    for (unsigned int iDac = 0; iDac < other->m_dacs.size(); iDac++) {
      m_dacs.push_back((other->m_dacs)[iDac]);
    }
  }

}

// $Header$
#ifndef CalibData_AncTaggerGain_h
#define CalibData_AncTaggerGain_h

#include "CalibData/RangeBase.h"

namespace CalibData {

  class AncTaggerGain : public RangeBase {
  public:
    AncTaggerGain(float val = 0.0, 
           unsigned isBad=0) : m_val(val), m_isBad(isBad) {}
    ~AncTaggerGain() {}

    float getVal() const {return m_val;}
    int getIsBad() const {return m_isBad;}
    virtual void update(RangeBase* other);

  private:
    float m_val;
    unsigned m_isBad;
  };
}

#endif

// $Header$
#ifndef CalibData_Gain_h
#define CalibData_Gain_h

#include "CalibData/Cal/RangeBase.h"

namespace CalibData {

  class Gain : public RangeBase {
  public:
    Gain(float gain = 0.0) : m_gain(gain) {}
    ~Gain() {}

    float getGain() {return m_gain;}

    virtual void update(RangeBase* other);

  private:
    float m_gain;

  };
}

#endif

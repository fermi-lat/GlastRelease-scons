// $Header$
#ifndef CalibData_AcdGain_h
#define CalibData_AcdGain_h

#include "CalibData/RangeBase.h"

namespace CalibData {

  class AcdGain : public RangeBase {
  public:
    AcdGain(float gain = 0.0, float sig = 0.0) : m_gain(gain), m_sig(sig) {}
    ~AcdGain() {}

    float getGain() const {return m_gain;}
    float getSig() const {return m_sig;}

    virtual void update(RangeBase* other);

  private:
    float m_gain;
    float m_sig;

  };
}

#endif

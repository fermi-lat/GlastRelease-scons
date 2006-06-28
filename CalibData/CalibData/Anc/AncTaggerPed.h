// $Header$
#ifndef CalibData_AncTaggerPed_h
#define CalibData_AncTaggerPed_h

#include "CalibData/RangeBase.h"

namespace CalibData {

  class AncTaggerPed : public RangeBase {
  public:
    AncTaggerPed(float val = 0.0, float rNoise = 0.0, float sNoise = 0.0,
           unsigned isBad=0) : m_val(val), m_rNoise(rNoise), 
                                m_sNoise(sNoise), m_isBad(isBad) {}
    ~AncTaggerPed() {}

    float getVal() const {return m_val;}
    float getRNoise() const {return m_rNoise;}
    float getSNoise() const {return m_sNoise;}
    int getIsBad() const {return m_isBad;}
    virtual void update(RangeBase* other);

  private:
    float m_val;
    float m_rNoise;
    float m_sNoise;
    unsigned m_isBad;
  };
}

#endif

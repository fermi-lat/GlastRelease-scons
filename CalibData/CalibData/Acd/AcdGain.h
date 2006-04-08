// $Header$
#ifndef CalibData_AcdGain_h
#define CalibData_AcdGain_h

#include "CalibData/RangeBase.h"

namespace CalibData {

  class AcdGain : public RangeBase {
  public:
    AcdGain(float peak = 0.0, float width = 0.0, unsigned status=0) :
      m_peak(peak), m_width(width), m_status(status) {}
    ~AcdGain() {}

    float getPeak() const {return m_peak;}

    float getWidth() const {return m_width;}
    unsigned getStatus() const {return m_status;}
    virtual void update(RangeBase* other);

  private:
    float m_peak;
    float m_width;
    unsigned m_status;
  };
}

#endif

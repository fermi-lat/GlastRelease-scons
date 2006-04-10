// $Header$
#ifndef CalibData_AcdPed_h
#define CalibData_AcdPed_h

#include "CalibData/RangeBase.h"

namespace CalibData {

  class AcdPed : public RangeBase {
  public:
    AcdPed(float mean = 0.0, float width = 0.0, unsigned status=0) :
      m_mean(mean), m_width(width), m_status(status) {}
    ~AcdPed() {}

    float getMean() const {return m_mean;}
    float getWidth() const {return m_width;}
    int getStatus() const {return m_status;}
    virtual void update(RangeBase* other);

  private:
    float m_mean;
    float m_width;
    int m_status;
  };
}

#endif

// $Header$
#ifndef CalibData_MuSlope_h
#define CalibData_MuSlope_h

#include "CalibData/RangeBase.h"

namespace CalibData {

  class MuSlope : public RangeBase {
  public:
    MuSlope(float slope = 0.0, float error = 0.0) : m_muSlope(slope),
    m_error(error) {}
    ~MuSlope() {}

    float getSlope() const {return m_muSlope;}
    float getError() const {return m_error;}

    virtual void update(RangeBase* other);

  private:
    float m_muSlope;
    float m_error;

  };
}

#endif

// $Header$
#ifndef CalibData_MuSlope_h
#define CalibData_MuSlope_h

#include "CalibData/Cal/RangeBase.h"

namespace CalibData {

  class MuSlope : public RangeBase {
  public:
    MuSlope(float slope = 0.0) : m_muSlope(slope) {}
    ~MuSlope() {}

    float getSlope() {return m_muSlope;}

    virtual void update(RangeBase* other);

  private:
    float m_muSlope;

  };
}

#endif

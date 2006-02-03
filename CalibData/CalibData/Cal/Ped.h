// $Header$
#ifndef CalibData_Ped_h
#define CalibData_Ped_h

#include "CalibData/RangeBase.h"

namespace CalibData {

  class Ped : public RangeBase {
  public:
    Ped(float av = 0.0, float sig = 0.0, float cosAngle=0.0) 
      : m_pedAvr(av), m_pedSig(sig), m_pedCosAngle(cosAngle) {}
    ~Ped() {}
    float getAvr() const {return m_pedAvr;}
    float getSig() const {return m_pedSig;}
    float getCosAngle() const {return m_pedCosAngle;}

    virtual void update(RangeBase* other);

  private:
    float m_pedAvr;
    float m_pedSig;
    float m_pedCosAngle;
  };
}
#endif

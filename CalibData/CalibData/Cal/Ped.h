// $Header$
#ifndef CalibData_Ped_h
#define CalibData_Ped_h

#include "CalibData/Cal/RangeBase.h"

namespace CalibData {

  class Ped : public RangeBase {
  public:
    Ped(float av = 0.0, float sig = 0.0) : m_pedAvr(av), m_pedSig(sig) {}
    ~Ped() {}
    float getAvr() {return m_pedAvr;}
    float getSig() {return m_pedSig;}

    virtual void update(RangeBase* other);

  private:
    float m_pedAvr;
    float m_pedSig;

  };
}
#endif

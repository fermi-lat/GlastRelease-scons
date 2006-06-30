// $Header$
#ifndef CalibData_AncQdcPed_h
#define CalibData_AncQdcPed_h

#include "CalibData/RangeBase.h"
#include <string>

namespace CalibData {

  class AncQdcPed : public RangeBase {
  public:
    AncQdcPed(float val = 0.0, float rms=0.0, unsigned isBad=0, 
              std::string dev="") : 
      m_val(val), m_rms(rms), m_isBad(isBad), m_device(dev) {}
    ~AncQdcPed() {}

    float getVal() const {return m_val;}
    float getRms() const {return m_rms;}
    int getIsBad() const {return m_isBad;}
    std::string getDevice() const {return m_device;}
    virtual void update(RangeBase* other);

  private:
    float m_val;
    float m_rms;
    unsigned m_isBad;
    std::string m_device;
  };
}

#endif

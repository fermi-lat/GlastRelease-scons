// $Header$
#ifndef CalibData_AcdVeto_h
#define CalibData_AcdVeto_h

#include "CalibData/RangeBase.h"

namespace CalibData {

  class AcdVeto : public RangeBase {
  public:
    AcdVeto(float veto = 0.0, float width = 0.0, unsigned status=0) :
      m_veto(veto), m_width(width), m_status(status) {}
    ~AcdVeto() {}

    float getVeto() const {return m_veto;}

    float getWidth() const {return m_width;}
    int getStatus() const {return m_status;}
    virtual void update(RangeBase* other);

  private:
    float m_veto;
    float m_width;
    int m_status;
  };
}

#endif

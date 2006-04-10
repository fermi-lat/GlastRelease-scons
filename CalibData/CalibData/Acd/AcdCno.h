// $Header$
#ifndef CalibData_AcdCno_h
#define CalibData_AcdCno_h

#include "CalibData/RangeBase.h"

namespace CalibData {

  class AcdCno : public RangeBase {
  public:
    AcdCno(float cno = 0.0, float width = 0.0, unsigned status=0) :
      m_cno(cno), m_width(width), m_status(status) {}
    ~AcdCno() {}

    float getCno() const {return m_cno;}

    float getWidth() const {return m_width;}
    int getStatus() const {return m_status;}
    virtual void update(RangeBase* other);

  private:
    float m_cno;
    float m_width;
    int m_status;
  };
}

#endif

// $Header$
#ifndef CalibData_LightAsym_h
#define CalibData_LightAsym_h

#include "CalibData/Cal/RangeBase.h"
#include <vector>

namespace CalibData {

  class LightAsym : public RangeBase {
  public:
    LightAsym(const std::vector<float>* values=0, float err = 0.0);
    ~LightAsym() {if (m_values) delete m_values;}

    const std::vector<float>* getValues() const {return m_values;}
    float getError() {return m_error;}

    virtual void update(RangeBase* other);

  private:
    std::vector<float>* m_values;
    float m_error;

  };
}

#endif

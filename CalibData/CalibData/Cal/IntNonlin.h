// $Header$
#ifndef CalibData_IntNonlin_h
#define CalibData_IntNonlin_h

#include "CalibData/Cal/RangeBase.h"
#include <vector>

namespace CalibData {

  /**  @class IntNonlin
    Class to keep track of integral non-linearity calibration data
    for a single range of a single face of a single crystal: just a
    set of values and a single error number.

    For each range (global to all crystals) there is a collection
    of dac values used to acquire data for that range; see
    @a DacCol.
  */
  class IntNonlin : public RangeBase {
  public:
    IntNonlin(const std::vector<float>* values=0, float err = 0.0);
    ~IntNonlin() {if (m_values) delete m_values;}

    const std::vector<float>* getValues() const {return m_values;}
    float getError() {return m_error;}

    virtual void update(RangeBase* other);

  private:
    std::vector<float>* m_values;
    float m_error;

  };
}

#endif

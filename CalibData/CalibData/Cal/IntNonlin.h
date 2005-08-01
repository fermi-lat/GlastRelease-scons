// $Header$
#ifndef CalibData_IntNonlin_h
#define CalibData_IntNonlin_h

#include "CalibData/RangeBase.h"
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
    IntNonlin(const std::vector<float>* values=0, float err = 0.0,
              const std::vector<float>* sdacs=0);
    ~IntNonlin() {
      if (m_values) delete m_values;
      if (m_sdacs) delete m_sdacs;
    }

    const std::vector<float>* getValues() const {return m_values;}
    const std::vector<float>* getSdacs() const {return m_sdacs;}
    float getError() {return m_error;}

    virtual void update(RangeBase* other);

  private:
    std::vector<float>* m_values;
    std::vector<float>* m_sdacs;
    float m_error;

  };
}

#endif

// $Header$
#ifndef CalibData_UniBase_h
#define CalibData_UniBase_h
#include "idents/TkrId.h"

namespace CalibData {
  using idents::TkrId;

  enum UniBy {
    UNI_BY_STRIP = 0,
    UNI_BY_GTFE = 1,
    UNI_BY_UNI  = 2
  };

  /** 
        Base class for per uniplane Tracker calibration data
   */
  class UniBase {

  public:
    UniBase(const TkrId& id=TkrId(), UniBy by=UNI_BY_STRIP, 
            unsigned nPer=64) 
      : m_by(by), m_stripPer(nPer)  {
      m_id.copy(id);
    }
    virtual ~UniBase() {}

    /// Derived classes will do a dynamic cast of argument, which 
    /// must be of same type, and then a deep copy.
    virtual void update(UniBase* ) {}

  protected:
    TkrId m_id;
    UniBy m_by;
    unsigned m_stripPer;  // only used when data is per gfte
  };

}
#endif

// $Header$

/** @class CalibBase
 *    Implementation of base class for all calibration data objects
 */

#include "CalibData/CalibBase.h"
#include "GaudiKernel/TimePoint.h" 
#include "GaudiKernel/StatusCode.h"

namespace CalibData {
  CalibBase::CalibBase() : m_validSince(0), m_validTill(0) {}

  CalibBase::CalibBase(const ITime& since, const ITime& till)
    m_validSince(0), m_validTill(0) {
    m_validSince = new TimePoint(since);
    m_validTill = new TimePoint(till);
}

  // Should be overridden by derived classes
  void CalibBase::update(CalibBase& obj) {
    delete m_validSince;
    delete m_validTill;

    till = new TimePoint(obj.validTill() );
    since = new TimePOint(obj.validSince() );
  }

  bool CalibBase::isValid() {
    return ((m_validSince != 0) && (m_validTill != 0)
            && (validSince() <= validTill())   );
  }

  bool Calibbase::isValid (const ITime& t) {
    if (!isValid()) return false;
    return validSince() <= t &&  t <= validTill();
  };

  const ITime& CalibBase::validSince() {
    return *m_validSince;
  }

  const ITime& CalibBase::validTill() {
    return *m_validTill;
  }

  void CalibBase::setValidity(const ITime& since, const ITime& till) {
    setValiditySince(since);
    setValidityTill(till);
  }

  void CalibBase::setValiditySince(const ITime& since) {
    delete m_validSince;
    m_validSince = new TimePoint(since);
  }

  void CalibBase::setValidityTill(const ITime& till) {
    delete m_validSince;
    m_validTill = new TimePoint(till);
  }

  StatusCode CalibBase::updateValidity() {
    return StatusCode::SUCCESS;
  }
}









}

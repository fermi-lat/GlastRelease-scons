// $Header$
#include "CalibData/CalibTime.h"

namespace {
  static longlong billion = 1000000000;
}
namespace CalibData {
  using facilities::Timestamp;

  CalibTime::CalibTime() : Timestamp() {}

  CalibTime::CalibTime(const Timestamp& stamp) :
    facilities::Timestamp(stamp) {}

  CalibTime::CalibTime(double julianDate) : facilities::Timestamp(julianDate)
  {}

  // This doesn't make any sense unless the input argument already "is"
  // of type CalibTime
  CalibTime::CalibTime(const ITime& time) {
    longlong absTime = time.absoluteTime();
    m_time = absTime / billion;
    m_nano = absTime - (m_time * billion);
  }

  ITime&  CalibTime::operator+=( const ITime& ) {
    /*   *this =  this->facilities::Timestamp::operator+=(CalibTime(other));
         return *this;
    */
    throw facilities::BadTimeInput
      ("CalibData::CalibTime Unsupported timestamp operation +=");
  }

  ITime&  CalibTime::operator-=( const ITime& ) {
    /*   *this =  this->facilities::Timestamp::operator-=(CalibTime(other));
         return *this;
    */
    throw facilities::BadTimeInput
      ("CalibData::CalibTime Unsupported timestamp operation -=");
  }


  ITime::AbsoluteTime CalibTime::absoluteTime() const {
    ITime::AbsoluteTime abs = m_time;
    abs *= billion;
    abs += m_nano;
    return abs;
  }

  ITime::DimensionedTime CalibTime::seconds() const {
    ITime::DimensionedTime dim = m_nano * 1000000000.0;
    dim += m_time;
    return dim;
  }

  bool CalibTime::operator==(const ITime& other) const {
    facilities::Timestamp me = *this;
    CalibTime otherCalib(other);
    return (me == otherCalib);
  }

  bool CalibTime::operator!=(const ITime& other) const {
    facilities::Timestamp me = *this;
    CalibTime otherCalib(other);
    return (me != otherCalib);
  }

  bool CalibTime::operator<=(const ITime& other) const {
    facilities::Timestamp me = *this;
    CalibTime otherCalib(other);
    return (me <= otherCalib);
  }

  bool CalibTime::operator>=(const ITime& other) const {
    facilities::Timestamp me = *this;
    CalibTime otherCalib(other);
    return (me >= otherCalib);
  }

  bool CalibTime::operator<(const ITime& other) const {
    facilities::Timestamp me = *this;
    CalibTime otherCalib(other);
    return (me < otherCalib);
  }

  bool CalibTime::operator>(const ITime& other) const {
    facilities::Timestamp me = *this;
    CalibTime otherCalib(other);
    return (me > otherCalib);
  }

  std::ostream& CalibTime::printOut(std::ostream& o) const {
    return o << getString();
  }
}

// $Header$

#ifndef CalibData_CalibTime_h
#define CalibData_CalibTime_h
/** @class CalibTime.h
    Has facilities::Timestamp as base and also keeps Gaudi::Time 
    representation internally.

    Note facilities::Timestamp, hence CalibTime, has get methods for 
    several different time formats, including julian date, ascii string, 
    and (seconds, nanoseconds) since Jan. 1 1970.
    Now (Oct. 2010) also adding Gaudi::Time as supported format
*/
#include "facilities/Timestamp.h"
#include "GaudiKernel/Time.h" 

namespace CalibData {
  class CalibTime : //public ITime,
                    public facilities::Timestamp
{

public:
  CalibTime();
  CalibTime(const facilities::Timestamp& stamp);
  CalibTime(double julianDate);
  
  // Need this one so that calibration objects can use CalibTime
  // objects internally to satisfy IValidity
  CalibTime(const Gaudi::Time& time);

  const Gaudi::Time& getGaudiTime() const;

  longlong seconds() const { return getGaudiTime().second(false); }
  int minutes() const { return getGaudiTime().minute(false); }
  int hours() const { return getGaudiTime().hour(false); }
  int days() const { return getGaudiTime().day(false); }

      
    std::ostream& printOut(std::ostream&) const;
private:
  Gaudi::Time m_gtime;

  };
}
#endif

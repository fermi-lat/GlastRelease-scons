// $Header$

#ifndef CalibData_CalibTime_h
#define CalibData_CalibTime_h
/** @class CalibTime.h
    Implements ITime interface and has facilities::Timestamp as base
    class.  Just need to implement a few more methods to satisfy
    ITime.
*/
#include "facilities/Timestamp.h"
#include "GaudiKernel/ITime.h" 

namespace CalibData {
  class CalibTime : public ITime,
                    public facilities::Timestamp
{

  public:
    CalibTime();
    CalibTime(const facilities::Timestamp& stamp);
    CalibTime(double julianDate);

    // Need this one so that calibration objects can use CalibTime
    // objects internally to satisfy IValidity
    CalibTime(const ITime& time);
    
    /// absoluteTime, seconds, etc. are needed to satisfy ITime interface
    ITime::AbsoluteTime absoluteTime() const;
    ITime::DimensionedTime seconds() const;
    ITime::DimensionedTime minutes() const {return seconds()/60.0;}
    ITime::DimensionedTime hours() const {return seconds()/3600.0;}
    ITime::DimensionedTime days() const {return seconds()/86400.0;}

    /// Reimplement operations for ITime interface using 
    /// facilities::Timestamp implementation
    bool              operator==(const ITime& other) const;

    bool              operator!=( const ITime& other) const;

    bool              operator<=( const ITime& other) const;

    bool              operator>=( const ITime& other) const;

    bool              operator< ( const ITime& other) const;

    bool              operator> ( const ITime& other) const;
  
    // adding   -- doesn't make sense for timestamps
    ITime&            operator+=(const ITime& other);
  
    // substraction  -- doesn't make sense for timestamps
    ITime&            operator-=(const ITime& other);
      
    std::ostream& printOut(std::ostream&) const;
  };
}
#endif

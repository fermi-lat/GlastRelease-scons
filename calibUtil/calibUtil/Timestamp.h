// $Header$
#ifndef CALIBUTIL_TIMESTAMP_H
#define CALIBUTIL_TIMESTAMP_H

#include <string>

namespace calibUtil {
  /** Encapsulation of timestamp for use as start/end of validity
      interval.  Resolution is a second.  Provide conversion to/from a 
      format that MySQL will accept for a datetime and provide comparison.

      Support the string format used for MySQL DATETIME type:

           yyyy-mm-dd hh:mm:ss

                          where 0 < mm < 13
                                0 < dd < 32
                               -1 < hh < 24
                               -1 < mm < 60 
                               -1 < ss < 60 

        o  only the first three fields are required.  Omitted trailing
           fields will be interpreted as equal to 0.
        o  by default : will be used to delimit fields, but user
           may specify an alternative in most circumstances
        o  leading zeros are optional

  */
  class Timestamp {
  public:
    /// Default constructor builds object representing current time
    Timestamp();

    /// Count seconds from the creation-according-to-unix, somewhere in 1970
    Timestamp(double seconds);

    /** Create a timestamp from an ascii string of standard form
                yyyy-mm-dd hh:mm:ss
        where only the first three fields are required.   

        If the string is invalid, object will represent unix creation time.

        Assume for now -- and preferably forever -- that all times are GMT.
    */
    Timestamp(const std::string& str);

    Timestamp(unsigned int year, unsigned int month, unsigned int day, 
              unsigned int hour = 0, unsigned int minute = 0, 
              unsigned int second = 0);

    std::string& timeString() {return m_strTime;}

    bool operator<(const Timestamp& other) {
      return m_time < other.m_time;
    }

    bool operator>(const Timestamp& other) {
      return m_time > other.m_time;
    }

    bool operator==(const Timestamp& other) {
      return m_time == other.m_time;
    }

    bool operator!=(const Timestamp& other) {
      return m_time != other.m_time;
    }

  private:

    static time_t toBinary(const std::string& strTime); 

    static void toString(time_t bin, std::string& strTime);

    /// internal binary rep of time
    time_t m_time;

    /// ascii equivalent
    std::string m_strTime;  
  };
}

#endif

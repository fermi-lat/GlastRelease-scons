// $Header$
#ifndef CALIBUTIL_TIMESTAMP_H
#define CALIBUTIL_TIMESTAMP_H

#include <string>

namespace calibUtil {
  /** Encapsulation of timestamp for use as start/end of validity
      interval.  Resolution is a second.  Provide conversion to/from a 
      format that MySQL will accept for a datetime and provide comparison.

      Support string format like this:

           yyyy:mm:dd:hh:mm:ss

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

    static const char DEFAULT_DELIM = ':';
    /** Create a timestamp from an ascii string with fields delimited by
        a one-character delimiter.  Fields, all represented with decimal
        strings, must start with year and decrease in significance.  Trailing
        fields after day may be omitted.  With default delimiter, valid
        input string look like 
               "2002:9:8"     (Sept. 8, 2002) or
               "2003:5:06:10" (May 6, 2003, 10 AM)
        If the string is invalid, object will represent unix creation time.

        Assume for now -- and preferably forever -- that all times are GMT.
    */
    Timestamp(const std::string& str, char delim = DEFAULT_DELIM);

    Timestamp(unsigned int year, unsigned int month, unsigned int day, 
              unsigned int hour = 0, unsigned int minute = 0, 
              unsigned int second = 0);

    std::string Timestamp::ascii(char delim = DEFAULT_DELIM);

    bool Timestamp::operator<(const Timestamp& other) {
      return m_time < other.m_time;
    }

    bool Timestamp::operator>(const Timestamp& other) {
      return m_time > other.m_time;
    }

    bool Timestamp::operator==(const Timestamp& other) {
      return m_time == other.m_time;
    }

    bool Timestamp::operator!=(const Timestamp& other) {
      return m_time != other.m_time;
    }

  private:

    static time_t toBinary(const std::string& strTime, 
                           char delim = DEFAULT_DELIM);

    static void toString(time_t bin, std::string& strTime, 
                        char delim = DEFAULT_DELIM);

    /// internal binary rep of time
    time_t m_time;

    /// For convenience also keep ascii rep using default delimiter
    std::string m_strTime;  
  };
}

#endif

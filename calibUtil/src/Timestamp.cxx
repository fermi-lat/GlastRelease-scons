// $Header$

#include <ctime>
#include <cstdlib>
#include <strstream>
#include "calibUtil/Timestamp.h"

namespace calibUtil {

  Timestamp::Timestamp(double seconds) : m_time(seconds) {
    toString(m_time, m_strTime);
  }

  Timestamp::Timestamp() {
    m_time = time(0);
    toString(m_time, m_strTime);
  }

  Timestamp::Timestamp(const std::string& str, char delim) {
    m_time = toBinary(str, delim);
    toString(m_time, m_strTime);
  }

  Timestamp::Timestamp(unsigned int year, unsigned int month, 
                       unsigned int day, unsigned int hour, 
                       unsigned int minute, unsigned int second) {
    struct tm fields;
    fields.tm_year = year - 1900;
    fields.tm_mon = month - 1;
    fields.tm_mday = day;
    fields.tm_hour = hour;
    fields.tm_min = minute;
    fields.tm_sec = second;
    fields.tm_wday = -1;
    fields.tm_yday = -1;
    
    fields.tm_isdst = 0;   // daylight savings time?

    m_time = timegm(&fields);
    toString(m_time, m_strTime);
  }

  time_t Timestamp::toBinary(const std::string& strTime, char delim) {
    // Break out fields
    struct tm fields;

    unsigned int pos;
    unsigned int oldPos = 0;

    pos = strTime.find(delim, oldPos);

    // Should have at least two occurrences of delim
    if (pos == strTime.size()) return 0;

    fields.tm_year = atoi((strTime.substr(oldPos, pos)).c_str()) - 1900;

    // on to month
    oldPos = pos + 1;
    pos = strTime.find(delim, oldPos);

    // Should have at least two occurrences of delim
    if (pos == strTime.size()) return 0;

    fields.tm_mon = atoi((strTime.substr(oldPos, pos)).c_str()) - 1;

    // day
    oldPos = pos + 1;

    pos = strTime.find(delim, oldPos);

    fields.tm_mday = atoi((strTime.substr(oldPos, pos)).c_str());

    // Remaining fields in string representation default to 0.
    fields.tm_hour = fields.tm_min = fields.tm_sec = 0;

    if (pos != strTime.size() ) {   // more fields to process
      oldPos = pos + 1;

      pos = strTime.find(delim, oldPos);

      fields.tm_hour = atoi((strTime.substr(oldPos, pos)).c_str());
      
      if (pos != strTime.size() ) {
        oldPos = pos + 1;
        pos = strTime.find(delim, oldPos);
        fields.tm_min = atoi((strTime.substr(oldPos, pos)).c_str());
        
        if (pos != strTime.size() ) {
          oldPos = pos + 1;
          pos = strTime.find(delim, oldPos);
          fields.tm_sec = atoi((strTime.substr(oldPos, pos)).c_str());
        }
      }
    }

    fields.tm_wday = -1;
    fields.tm_yday = -1;
    fields.tm_isdst = 0;
    return timegm(&fields);
  }
  
  void Timestamp::toString(time_t bin, std::string& strTime, 
                           char delim) {
    struct tm * fields = gmtime(&bin);

    std::string dstring(&delim, 1);
    strTime.resize(0);
    std::strstream s;
    s << fields->tm_year << dstring << fields->tm_mon << dstring;
    s << fields->tm_mday << dstring << fields->tm_hour << dstring;
    s << fields->tm_min << dstring << fields->tm_sec;
    strTime += s.str();

  }


}

// $Header$

#include <ctime>
#include <cstdlib>
#include <strstream>
#include "calibUtil/Timestamp.h"

namespace calibUtil {

  Timestamp::Timestamp(double seconds) : m_time((time_t) seconds) {
    toString(m_time, m_strTime);
  }

  Timestamp::Timestamp() {
    m_time = time(0);
    toString(m_time, m_strTime);
  }

  Timestamp::Timestamp(const std::string& str) {
    m_time = toBinary(str);
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

  time_t Timestamp::toBinary(const std::string& strTime) {
    // Break out fields
    struct tm fields;

    unsigned int pos;
    unsigned int oldPos = 0;

    // Three delimiter characters may occur in ascii time representation.
    // First is hyphen.
    char  delim = '-';

    pos = strTime.find(delim, oldPos);

    // Must have two occurrences of hyphen delimiter
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

    // A space separate time from date (if time is present at all)
    delim = ' ';
    pos = strTime.find(delim, oldPos);

    fields.tm_mday = atoi((strTime.substr(oldPos, pos)).c_str());

    // Remaining fields in string representation default to 0.
    fields.tm_hour = fields.tm_min = fields.tm_sec = 0;

    if (pos != strTime.size() ) {   // more fields to process
      delim = ':';
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
  
  void Timestamp::toString(time_t bin, std::string& strTime) {
    struct tm * fields = gmtime(&bin);

    strTime.resize(0);
    std::strstream s;
    s << (fields->tm_year + 1900) << "-" << (fields->tm_mon + 1) << "-";
    s << fields->tm_mday << " " << fields->tm_hour << ":";
    s << fields->tm_min << ":" << fields->tm_sec;
    strTime += s.str();

  }


}

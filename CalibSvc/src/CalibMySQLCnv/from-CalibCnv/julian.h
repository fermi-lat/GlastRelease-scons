#ifndef JULIANDATE_H
#define JULIANDATE_H


/** 
     Convert dates starting with Jan 1 2000 and thereafter between
     Julian date rep.  and (year, month, day, hour, minute, sec)
     Time part of the representation is gmt.

     Components of time are always 0-based; that is, months range
     from 0 to 11, hours from 0 to 23, etc.
*/ 
class juliandate {

public:
  juliandate(double jul = JULAN2000);
  juliandate(int year, int month, int day, int hour=0, int minute=0,
             double sec = 0.0);

  double getJulian() {return m_julian;}

  int getYear() {return m_year;}
 
  /// Return zero-based month
  int getMonth() {return m_month;}

  /// Return zero-based day of month
  int getDay() {return m_day;}

  /// Return zero-based hour (gmt)
  int getHours() {return m_hour;}

  int getMinutes() {return m_minute;}

  double getSeconds() {return m_sec;}


private:
  /// Julian date for start of Jan 1, 2000
  static const double julian2000 = 2451544.5;
  static const int    secPerDay =  (24*60*60);

  double m_julian;
  int    m_year;
  int    m_month;
  int    m_day;
  int    m_hour;
  int    m_minute;
  double m_sec;
}
#endif

// $Header$
#ifndef EVENT_TIMESTAMP_H
#define EVENT_TIMESTAMP_H 1


// Include files
#include <iostream>
#include "GaudiKernel/StreamBuffer.h"
#include "Event/TopLevel/Definitions.h"

/*! \class TimeStamp
\brief encapsulate the time.
Elapsed? absolute? Currently a double, in units of seconds.
one event. 
*/

class dateRep{
public:
    dateRep():m_years(0),m_months(0),m_days(0),m_hours(0){}
    
    dateRep(int years, int months, int days, double hours):
    m_years(years),m_months(months),m_days(days),m_hours(hours)
    {}
    int years() const{return m_years;}
    int months() const{return m_months;}
    int days() const{return m_days;}
    double hours() const{return m_hours;}
    void reSet(int years, int months, int days, double hours){
        m_years = years;
        m_months = months;
        m_days = days;
        m_hours = hours;
    }
    int m_years;
    int m_months;
    int m_days;
    double m_hours;
};

class TimeStamp                                                                {
    
public:
    
    /// Constructors
    TimeStamp()
        : m_time(0), m_startSimTime(GetJD(2006.,7,18,0.0)),m_startSimDate(2006.,7,18,0.0)
    { GetCurrentDate();}
    
    TimeStamp( double t )
        : m_time(t), m_startSimTime(GetJD(2006.,7,18,0.0)),m_startSimDate(2006.,7,18,0.0)
    { GetCurrentDate();}
    
    /// initialize with a date
    TimeStamp(int years,int months,int days,double hours)
        :m_startSimTime(GetJD(2006.,7,18,0.0)),m_startSimDate(2006.,7,18,0.0),
        m_time(GetJD(years,months,days,hours) - m_startSimTime){
        m_time = GetJD(years,months,days,hours) - m_startSimTime;
        GetCurrentDate();}
    
    ///initialize with a current date and a start date.
    TimeStamp(int years,int months,int days,double hours,int yearsBegin, int monthsBegin, int daysBegin, double hoursBegin)
        :m_startSimTime(GetJD(yearsBegin,monthsBegin,daysBegin,hoursBegin)),m_startSimDate(yearsBegin,monthsBegin,daysBegin,hoursBegin),
        m_time(GetJD(years,months,days,hours) - m_startSimTime){
        m_time = GetJD(years,months,days,hours) - m_startSimTime;
        GetCurrentDate();}
    
    ///Destructor
    ~TimeStamp(){}
    
    /// Retrieve time
    double time() const{
        return m_time;
    }
    
    ///retreive date
    dateRep date() const{
        return m_date;
    }
    
    //Retreive the starting date
    dateRep startDate() const{
        return m_startSimDate;
    }
    
    /// Update time 
    void setTime( double value) {
        m_time = value;
        GetCurrentDate();
    }
    
    void setTime(int years,int months,int days,double hours){
        m_time = GetJD(years,months,days,hours);
        GetCurrentDate();
    }
    
    operator double()const { return time(); }
    
    /// Serialize the object for writing
    friend StreamBuffer& operator<< ( StreamBuffer& s, const TimeStamp& obj )    {
        return s << obj.m_time;
    }
    /// Serialize the object for reading
    friend StreamBuffer& operator>> ( StreamBuffer& s, TimeStamp& obj )          {
        return s >> obj.m_time;
    }
    
    /// Output operator (ASCII)
    friend std::ostream& operator<< ( std::ostream& s, const TimeStamp& obj )    {
        return obj.fillStream(s);
    }
    /// Fill the output stream (ASCII)
    std::ostream& fillStream( std::ostream& s ) const{
        return s << "class TimeStamp : "
            << " Year: " << m_date.m_years
            << " Month: " << m_date.m_months
            << " Day: " << m_date.m_days
            << " Hour: " << m_date.m_hours << " , "
            << EventField( EventFormat::field12 )
            << m_time;
    }
    
private:
    
    void GetCurrentDate(){
        dateRep returned;
        
        int years,months,days;
        double hours, secondsLeft;
        secondsLeft = m_time;
        years = m_startSimDate.years();
        months = m_startSimDate.months();
        days = m_startSimDate.days();
        hours = m_startSimDate.hours();
        double secondsPerMonth = 86400*30.4375/*6001*/;
        double secondsPerYear = 86400*365.25;
        double secondsPerDay = 86400;
        double secondsPerHour = 3600;
        
        
        //[SET DATE VALUES FROM TIME HERE]
        
        while(secondsLeft >= secondsPerYear){
            years++;
            secondsLeft -= secondsPerYear;
        }
        while(secondsLeft >= secondsPerMonth){
            months++;
            secondsLeft -= secondsPerMonth;
            if(months > 12){
                months -=12;
                years++;}
        }
        while(secondsLeft >= secondsPerDay){
            days++;
            secondsLeft -= secondsPerDay;
            if(days > 30){
                days -=30;
                months++;
                if(months > 12){
                    months -=12;
                    years++;
                }
            }
        }
        while(secondsLeft >= secondsPerHour){
            
            hours++;
            secondsLeft -= secondsPerHour;
        }
        hours += secondsLeft/secondsPerHour;
        if(hours >=24){
            hours -=24;
            days++;
            if(days > 30){
                days -=30;
                months++;
                if(months > 12){
                    months -=12;
                    years++;
                }
            }
        }
            
            m_date.reSet(years,months,days,hours);
            //return returned;
        }
        
        double GetJD(int years,int months,int days,double hours)
        {
            double A,B,D;
            double J_D;
            double secondsPerDay = 86400;
            /*long int*/double C;
            
            if (months > 2);
            else {
                years = years - 1;
                months = months + 12;
            }
            A = (years / 100); B = 2 - A + (A / 4);
            C = /*(long int)*/(365.25 * years); if (years < 0) C = C - 1;
            D = /*(int)*/(30.4375/*6001*/ * (months + 1));
            J_D = B + C + D + days + 1720994.5+ hours / 24.;
            return J_D*secondsPerDay;
        }   /* daysrno_Giuliano */    
        
        /// Time
        double m_time; //this is the offset number of seconds since the simulation began
        dateRep m_date;
        double m_startSimTime;// this is the number corresponding to the time of the simulation beginning.
        dateRep m_startSimDate;
    };
    
    
#endif    // LHCBEVENT_TIMESTAMP_H

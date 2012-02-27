/********************************************************************************
* Author:  Christoph Singewald						             				*
* develop@nuvitec.com								             				*
* www.nuvitec.com																*
* see FXDate.h										             				*
*********************************************************************************/

#include "FXDate.h"
 

namespace FX {
	
	static int	DaysSoFar[][13] =
	{
		{0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365},
		{0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}
	};
	
	static int MONTH_LENGTH[] = {31,28,31,30,31,30,31,31,30,31,30,31}; // 0-based
	static int LEAP_MONTH_LENGTH[] = {31,29,31,30,31,30,31,31,30,31,30,31}; // 0-based
	
	/* default language is english 
	/  temp hack, not the best way :)
	/  when FOX has international support we can change that.
	*/
 
#ifdef _GERMAN
	// German names
	static const char* SHORT_DAY_NAMES[] = {"SO","MO","DI","MI","DO","FR","SA"}; // 0-based
	static const char* MONTH_NAMES[]
		= {"Januar","Februar","März","April","Mai","Juni","Juli","August","September","Oktober","November","Dezember"}; // 0-based
	static const char* DAY_NAMES[]
		= {"Sonntag","Montag","Dienstag","Mittwoch","Donnerstag","Freitag","Samstag"}; // 0-based
	static const char *quater = "Quartal";
#else 
#ifdef _FRENCH // by Marc Duren marc.duren@9online.fr
    static const char* MONTH_NAMES[]=  {"Janvier","Février","Mars","Avril","Mai","Juin","Juillet","Aout","Septembre","Octobre","Novembre","Décembre"}; // 0-based
	static const char* DAY_NAMES[] =  {"Dimanche","Lundi","Mardi","Mercredi","Jeudi","Vendredi","Samedi"}; // 0-based
	static const char *quater = "Trimestre";
	static const char* SHORT_DAY_NAMES[] =  {"DIM","LUN","MAR","MER","JEU","VEN","SAM"}; // 0-based

#else
	// English 
	static const char* SHORT_DAY_NAMES[] = {"SUN","MON","TUE","WED","THU","FRI","SAT"}; // 0-based
	static const char* MONTH_NAMES[]
		= {"January","February","March","April","May","June","July","August","September","October","November","December"}; // 0-based
	static const char* DAY_NAMES[]
		= {"Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"}; // 0-based
	static const char *quater = "quater";
#endif
#endif	
	
	/* Constructors */
	
	FXDate::FXDate() 
	{
		
		tm			*local_time;
		time_t timer	   = time(NULL);
		local_time = localtime(&timer);
		
		month = (short)(local_time->tm_mon + 1);
		day   = (short) local_time->tm_mday;
		year  = (short)(local_time->tm_year + 1900);
		hour = (short)(local_time->tm_hour);
		minute = (short)(local_time->tm_min);
		second= (short)(local_time->tm_sec);
		
		mdy_to_julian ();
	};
	
	FXDate::FXDate (const FXDate &dt)
	{
		month = dt.month;
		day   = dt.day;
		year  = dt.year;
		hour = dt.hour;
		minute = dt.minute;
		second= dt.second;
		
		mdy_to_julian();
	}
	
	FXDate::FXDate(const time_t time) 
	{ 
		
		tm	*local_time;
		
		local_time = localtime(&time);
		
		month = (short)(local_time->tm_mon + 1);
		day   = (short) local_time->tm_mday;
		year  = (short)(local_time->tm_year + 1900);
		hour = (short)(local_time->tm_hour);
		minute = (short)(local_time->tm_min);
		second= (short)(local_time->tm_sec);
		
		mdy_to_julian ();
	};
	
	
	FXDate::FXDate(int year, int month, int day, int hour, int minute, int second)
	{
		setNew(year, month, day, hour, minute,second);
	}
	
	FXDate::FXDate(int quarter)
	{
		setQuarter(quarter);
	}
	FXDate::FXDate(FXString datestring)
	{
		setDate(datestring);
	}
	/***
	IsLaepYear()
	
	 http://www.kalenderlexikon.de/anzeigen.php?Eintrag=CalcSchaltjahr
	 
	  Im Gregorianischen Kalender (ab dem 15.Oktober 1582) ist die Definition folgendermaßen:
	  Alle 4 Jahre ein Schaltjahr, alle vollen 100 Jahre nicht,
	  volle 400 Jahre jedoch wieder. Dies bedeutet,
	  der Gregorianische Kalender hat 97 Schaltjahre in 400 Jahren. Das ergibt in 3300 Jahren einen Fehler
	  von einem Tag im Vergleich zum tropischen Jahr.
	  
	****/
	bool FXDate::isLeapYear()
	{
		
		julian_to_mdy();
		int y = year;
		
		if(y%4==0 && y%100==0)
			return (y%400==0);
		else
			return (y%4==0);
		
		
	};
	
	
	bool FXDate::checkStringFormat(FXString datestring)
	{
		
		if(datestring.length()>10)
			return false;
		
		unsigned int day=0, month=0, year=0;
		
		if (3 != datestring.scan("%u.%u.%u", &day, &month, &year))
			return false;
		
		if(day > 31 || month>12 || year<999 || year>9999)
			return false;
		
		
		
		
		return true;
		
	}
	
	
	/* Set the Elements */
	void FXDate::setDate(FXString datestring)
	{
		
		if(datestring.length()>10)
			return;
		
		unsigned int day=0, month=0, year=0;
		
		if (3 != datestring.scan("%u.%u.%u", &day, &month, &year))
			return;
		
		if(month>12 || year<999 || year>9999)
			return;
		
		
		setNew(year,month,day);
		
	}
	
	void FXDate::setNew(int YY, int MM, int DD, int hh, int mm, int ss)
	{
		
		month=MM;
		day=DD;
		year=YY;
		FXASSERT(month > 0);
		FXASSERT(month < 13);
		FXASSERT(day   > 0);
		FXASSERT(day   < 32);
		
		hour=hh;
		minute=mm;
		second=ss;
		mdy_to_julian();	
	}
	
	
	void FXDate::setSystime(time_t systime)
	{
		tm*  local_time = localtime(&systime);
		
		if(!local_time) return;
		month = (short)(local_time->tm_mon + 1);
		day   = (short) local_time->tm_mday;
		year  = (short)(local_time->tm_year + 1900);
		hour = (short)(local_time->tm_hour);
		minute = (short)(local_time->tm_min);
		second= (short)(local_time->tm_sec);
		
		mdy_to_julian ();
		
	}
	
	void FXDate::setYear(int YYYY)
	{
		
		year=YYYY;
		mdy_to_julian ();
		
		
	}
	
	void FXDate::setMonth(int MM)
	{ 
		month=MM;
		if(month<1)
		{
			month=12;
			year--;
		}
		else if(month>12)
		{
			month=1;
			year++;
		}

		mdy_to_julian ();
		
		
	}
	
	
	void FXDate::setDay(int DD)
	{
		day=DD;
		mdy_to_julian ();
		
	}
	
	void FXDate::setHour(int hh)
	{
		
		hour=hh;
		mdy_to_julian ();
		
		
		
	}
	void FXDate::setMinute(int mm)
	{
		
		minute=mm;
		mdy_to_julian ();
		
		
		
	}
	void FXDate::setSec(int ss)
	{
		
		second=ss;
		mdy_to_julian ();
		
		
		
	}
	
	void FXDate::setHMS(int hh, int mm, int ss)
	{
		second=ss;
		minute=mm;
		hour=hh;
		mdy_to_julian ();
		
	}
	
	void FXDate::setQuarter(int quarter)
	{
		
		int dy=(quarter<1)?(quarter/4)-1:(quarter-1)/4;
		int yy = year+dy;
		
		quarter = (quarter-1)%4;
		quarter=(quarter<0)?quarter+=5:quarter+1;
		
		switch(quarter)
		{
		case 1:
			setNew(yy,01,01,1,0);
			break;
		case 2:
			setNew(yy,04,01,1,0);
			
			break;
		case 3:
			setNew(yy,07,01,1,0);
			break;
		case 4:
			setNew(yy,10,01,1,0);
			break;
		}
		
		mdy_to_julian();
	}
	
	time_t FXDate::getSystime() 
	{ 
		struct tm tmRep;
		time_t t;
		
		julian_to_mdy();
		
		int y=year;
		int m=month;
		int d=day;;
		
		if ( y < 1970 )
		{
			y = 70;
			m = 1;
			d = 1;
		}
		tmRep.tm_year = y - 1900 ;
		tmRep.tm_mon = m-1;
		tmRep.tm_mday = d;
		tmRep.tm_hour = hour;
		tmRep.tm_min = minute;
		tmRep.tm_sec = second;
		tmRep.tm_isdst = 0;
		
		t = mktime( &tmRep );
		return t;
		
	};
	
	
	int FXDate::getYear()
	{
		julian_to_mdy();
		return year;
	};
	int FXDate::getMonth()
	{
		julian_to_mdy();
		return month;
		
	};
	int FXDate::getDay()
	{
		julian_to_mdy();
		return day;
	};
	
	int FXDate::getHour()
	{
		julian_to_mdy();
		return hour;
		
	};
	int FXDate::getMinute()
	{
		julian_to_mdy();
		return minute;
		
	};
	int FXDate::getSecond()
	{
		julian_to_mdy();
		return second;
	};
	int FXDate::getDayOfWeek() {
		
		julian_to_mdy();
		return day_of_week-1;
		
	};
	
	int FXDate::getDaysOfMonth(int month)
	{
		if(isLeapYear())
			return LEAP_MONTH_LENGTH[month-1];
		else
			return MONTH_LENGTH[month-1];
	};
	
	
	FXString FXDate::getDayName()
	{
		
		return DAY_NAMES[this->getDayOfWeek()];
	}
	FXString FXDate::getShortDayName()
	{
		
		return SHORT_DAY_NAMES[this->getDayOfWeek()];
	}
	int FXDate::getTTMM()
	{
		int ttmmjj = 0;
		
		julian_to_mdy();
		
		ttmmjj = day* 100;
		ttmmjj += month;
		
		
		return ttmmjj;
		
	}
	int FXDate::getTTMMJJ()
	{
		int ttmmjj = 0;
		
		julian_to_mdy();
		
		ttmmjj = day * 10000;
		ttmmjj += month * 100;
		
		
		ttmmjj += year - (abs(year/100) * 100);
		
		return ttmmjj;
		
	}
	
	short FXDate::getYY()
	{
		julian_to_mdy();
		return year - (abs(year/100) * 100);
	}
	
	// Operators
	
	void  FXDate::operator= (FXDate  const &param) {
		this->julian  = param.julian;
		julian_to_mdy();
	}
	
	
	FXDate FXDate::operator- (FXDate param) {
		FXDate temp;
		temp.julian =julian-param.julian;;
		julian_to_mdy();
		return (temp);
	}
	FXDate FXDate::operator++(int)  {
		day++;
		if(day<getDaysOfMonth(month)) {
			month++;
			if(month<12)
				year++;
			day=1;
		}
		mdy_to_julian();
		return *this;
	}

	FXDate FXDate::operator--(int)  {
		day--;
		if(day<1) {
			month--;
			if(month<1)
				year--;
			day=getDaysOfMonth(month);
		}
		mdy_to_julian();
		return *this;
	}


	FXDate FXDate::operator+ (FXDate param) {
		FXDate temp;
		temp.julian =julian+param.julian;;
		julian_to_mdy();
		return (temp);
	}
	
	bool FXDate::operator > (FXDate param)
	{
		
		long  dif =param.julian-julian;
		return (dif<0)?true:false;
	}
	bool FXDate::operator < (FXDate param)
	{
		long  dif =param.julian-julian;
		return (dif>0)?true:false;
	}
	bool FXDate::operator >= (FXDate param)
	{
		
		long  dif =param.julian-julian;
		return (dif<=0)?true:false;
	}
	bool FXDate::operator <= (FXDate param)
	{
		long  dif =param.julian-julian;
		return (dif>=0)?true:false;
	}
	
	bool FXDate::operator == (FXDate param)
	{
		if(year==param.year &&
			month==param.month && 
			day==param.day)
			return true;
		else
			return false;
	}
	
	int FXDate::getQuarter()
	{
		
		julian_to_mdy();
		
		FXDate begin_q1(year,01,01,0,0);
		FXDate begin_q2(year,04,01,0,0);
		FXDate begin_q3(year,07,01,0,0);
		FXDate begin_q4(year,10,01,0,0);
		
		
		if(*this>=begin_q1 && *this<begin_q2) {
			return 1;
		}
		else
			if(*this>=begin_q2 && *this<begin_q3) {
				return 2;
			}
			else
				if(*this>=begin_q3 && *this <begin_q4) {
					return 3;
				}
				else
					if(*this>=begin_q4) {
						return 4;
						
					}
					
					return -1;
					
	}
	FXString FXDate::text()
	{
		FXString buffer;
		buffer.format("%02d.%02d.%02d",day,month,year);
		return buffer;
	}
	FXString FXDate::quarterToShortText()
	{
		
		julian_to_mdy();
		
		FXDate begin_q1(year,01,01,0,0);
		FXDate begin_q2(year,04,01,0,0);
		FXDate begin_q3(year,07,01,0,0);
		FXDate begin_q4(year,10,01,0,0);
		
		
		FXString buffer;
		
		if(*this>=begin_q1 && *this<begin_q2) {
			buffer.format("Q1/%d",year);
		}
		else
			if(*this>=begin_q2 && *this<begin_q3) {
				buffer.format("Q2/%d",year);
			}
			else
				if(*this>=begin_q3 && *this <begin_q4) {
					buffer.format("Q3/%d",year);
				}
				else
					if(*this>=begin_q4) {
						buffer.format("Q4/%d",year);
						
					}
					
					
					return buffer;
					
	}
	FXString FXDate::quarterToText()
	{
		
		
		julian_to_mdy();
		
		FXDate begin_q1(year,01,01,0,0);
		FXDate begin_q2(year,04,01,0,0);
		FXDate begin_q3(year,07,01,0,0);
		FXDate begin_q4(year,10,01,0,0);
		
		
		FXString buffer;
		
		if(*this>=begin_q1 && *this<begin_q2) {
			buffer.format("1.%s: 01.02.%d - 31.03.%d",quater,year,year);
		}
		else
			if(*this>=begin_q2 && *this<begin_q3) {
				buffer.format("2.%s: 01.04.%d - 30.06.%d",quater,year,year);
			}
			else
				if(*this>=begin_q3 && *this <begin_q4) {
					buffer.format("3.%s: 01.07.%d - 30.09.%d",quater,year,year);
				}
				else
					if(*this>=begin_q4) {
						buffer.format("4.%s: 01.10.%d - 31.12.%d",quater,year,year);
						
					}
					
					return buffer;
					
	}
	
	FXString  FXDate::textWithDayName()
	{
		FXString buffer;
		buffer.format("%s  %02d.%02d.%04d",getDayName().text(),day,month,year);
		return buffer;
	}
	
	
	
	FXString  FXDate::textWidthoutYear()
	{
		
		FXString buffer;
		buffer.format("%s  %02d.%02d",getDayName().text(),day,month);
		return buffer;
	}
	
	FXString  FXDate::timeToText()
	{
		
		FXString buffer;
		buffer.format("%02d:%02d",hour,minute);
		return buffer;   
		
	}
	FXString FXDate::monthYearToText()
	{
		FXString buffer;
		buffer.format("%s %02d",MONTH_NAMES[month-1],year);
		return buffer;
	}
	
	FXString FXDate::getMonthname()
	{
		return MONTH_NAMES[month-1];
	}
	FXString FXDate::getShortDayName(FXint index)
	{
		if(index>6) return DAY_NAMES[0];
		return SHORT_DAY_NAMES[index];
	}
	
	
	FXString FXDate::getDayName(FXint index)
	{
		if(index>6) return DAY_NAMES[0];
		return DAY_NAMES[index];
	}
	
	FXString FXDate::getMonthName(FXint index)
	{
		if(index>11) return MONTH_NAMES[11];
		return MONTH_NAMES[index];
	}
	
	
	/*************************************************************
	* Conversion routines
	**************************************************************/
	
	void FXDate::julian_to_wday (void)
	{
		day_of_week = (enum Wday) ((julian + 2) % 7 + 1);
	}
	
	/*************************************************************/
	
	
#define OCT5_1582		(2299160L)		/* "really" 15-Oct-1582  */
#define OCT14_1582		(2299169L)		/* "really"  4-Oct-1582  */
#define JAN1_1			(1721423L)
	
#define YEAR			(365)
#define FOUR_YEARS		(1461)
#define CENTURY 		(36524L)
#define FOUR_CENTURIES	(146097L)
	
	void FXDate::julian_to_mdy ()
	{
		long	z,y;
		short 	m,d;
		int 	lp;
		
		z = julian+1;
		if (z >= OCT5_1582)
		{
			z -= JAN1_1;
			z  = z + (z/CENTURY)  - (z/FOUR_CENTURIES) -2;
			z += JAN1_1;
			
		}
		
		z = z - ((z-YEAR) / FOUR_YEARS);		/* Remove leap years before current year */
		y = z / YEAR;
		
		d = (short) (z - (y * YEAR));
		
		y = y - 4712;				/* our base year in 4713BC */
		if (y < 1)
			y--;
		
		lp = !(y & 3);				/* lp = 1 if this is a leap year. */
		
		if (d==0)
		{
			y--;
			d = (short) (YEAR + lp);
		}
		
		m  = (short) (d/30);		/* guess at month */
									
		while (DaysSoFar[lp][m] >=d)
			m--;					/* Correct guess. */
		
		d = (short) (d - DaysSoFar[lp][m]);
		
		day = d;
		
		month = (short) (m+1);
		
		year = (short) y;
		
		julian_to_wday ();
	}
	
	/*************************************************************
	
	 The original here was far more complicated then it needed to be.
		 What we need to keep in mind is the simple rule:
		  Before 10/4/1585, a leap year occured every 4 years.
		  After  10/15/1585, leap years were skipped on centuries
			not divisible by 400. Plus 10 days were skipped to adjust
			for the past error.
	
   	*************************************************************/

	void FXDate::mdy_to_julian (void)
	{
		int		a;
		int		work_year=year;
		long	j;
		int 	lp;
		
		/* correct for negative year  (-1 = 1BC = year 0) */
		
		if (work_year < 0)
			work_year++;
		
		lp = !(work_year & 3);			/* lp = 1 if this is a leap year. */
		
		j =
			((work_year-1) / 4)		+		/* Plus ALL leap years */
			DaysSoFar[lp][month-1]	+
			day					+
			(work_year * 365L)	+		/* Days in years */
			JAN1_1 			+
			-366;						/* adjustments */
		
		/* deal with Gregorian calendar */
		if (j >= OCT14_1582)
		{
			
			a = (int)(work_year/100);
			j = j+ 2 - a + a/4;			/* Skip days which didn't exist. */
		}
		
		julian = j;
		
		julian_to_wday ();
	}

	void FXDate::validate()
	{
		
	}
	
}

/********************************************************************************
*                                                                               *
*             		F X D a t e    W i d g e t	                        		*
*                                                                               *
*********************************************************************************
* Author  Christoph Singewald						             				*
* develop@nuvitec.com								             				*
* www.nuvitec.com									             				*
*********************************************************************************
*********************************************************************************
* The FoxToolkit library is														*
* Copyright (C) 1998,2003 by Jeroen van der Zijp.   All Rights Reserved.        *
*********************************************************************************
* This library is free software; you can redistribute it and/or                 *
* modify it under the terms of the GNU Lesser General Public                    *
* License as published by the Free Software Foundation; either                  *
* version 2.1 of the License, or (at your option) any later version.            *
*                                                                               *
* This library is distributed in the hope that it will be useful,               *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU             *
* Lesser General Public License for more details.                               *
*                                                                               *
* You should have received a copy of the GNU Lesser General Public              *
* License along with this library; if not, write to the Free Software           *
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.    *
*********************************************************************************
* $Id$   	                    *
*********************************************************************************
*********************************************************************************
* FXDate is	based on source of James M. Curran									*
* 																				 *
*********************************************************************************/
/*
  Notes: Please report errors to: develop@nvutiec.com
		
*/
/*

  ToDo:
	
	 * set and get methodes with formatstrings like printf,scanf ...
     * bether documentation
	 * testing and testing ....
	 * better international support


*/


#ifndef _FXDATE_H
#define _FXDATE_H


#include <ctime>
#include "fx.h"

namespace FX { 
	
	class FXDate {
		
	private:
		
		unsigned long julian;
		short 	year;		 
		short	month;			 
		short 	day;			 
		char 	day_of_week;	/* 	1 = Sunday, ... 7 = Saturday */
		char	separator;
		short   hour;
		short   minute;
		short   second;
		
		FXint lang;
		
	private:
		void julian_to_mdy ();         
		void julian_to_wday ();        
		void mdy_to_julian ();   
		
		void validate();
		
		
	public:
		
		enum Wday {SUN=1,MON,TUE,WED,THU,FRI,SAT};
		
		FXDate();
		FXDate(const FXDate &dt);
		FXDate(const time_t time);
		FXDate(int year, int month, int day,int hour=0, int minute=0, int second=0);
		FXDate(FXString datestring);
		FXDate(int quarter);
		
		~FXDate() {}
		void operator= (FXDate const  &param);
		FXDate operator + (FXDate param);
		FXDate operator - (FXDate param);
		FXDate operator ++(int);
		FXDate operator --(int);
		bool operator < (FXDate param);
		bool operator > (FXDate param);
		bool operator <= (FXDate param);
		bool operator >= (FXDate param);
		bool operator == (FXDate param);
		bool operator != (FXDate param) { return !(*this==param);};
		
		inline  unsigned long juliandate() { return julian; };
		
		time_t getSystime();
		int getYear();
		int getMonth();
		int getDay();
		int getDayOfWeek();			/* 	0 = Sunday, ... 6 = Saturday */
		int getDaysOfMonth(int month);
		int getHour();
		int getMinute();
		int getSecond();
		FXString getDayName();
		FXString getShortDayName();
		int getTTMM();
		int getTTMMJJ();
		short getYY();
		int getQuarter();
		
		void setYear(int YYYY);
		void setMonth(int MM);
		void setDay(int DD);
		void setMinute(int mm);
		void setHour(int hh);
		void setSec(int ss);
		void setNew(int YY, int MM, int DD,int hh = 0, int mm=0, int ss = 0);
		void setSystime(time_t time);
		void setDate(FXString newdate);
		void setQuarter(int quarter);
		void setHMS(int hh, int mm, int ss);
		
		bool isLeapYear();
		
		FXString text();
		FXString textWithDayName();
		FXString textWidthoutYear();
		FXString timeToText();
		FXString monthYearToText();
		FXString quarterToShortText();
		FXString quarterToText();
		FXString getMonthname();
		
		
		
		static bool checkStringFormat(FXString datestring);
		
		static FXString getShortDayName(FXint index); // 0 ... Mon 
		static FXString getDayName(FXint index); // 0 ... Monday 
		static FXString getMonthName(FXint index); // 0 ... January
		
		
		
		
	};
}


#endif


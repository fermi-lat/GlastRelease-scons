/********************************************************************************
*                                                                               *
*                       F X D a t e P i c k e r   W i d g e t                   *
*                                                                               *
*********************************************************************************
* Author:  Christoph Singewald                                                  *
* develop@nuvitec.com                                                           *
* www.nuvitec.com                                                               *
*********************************************************************************
*********************************************************************************
* The FoxToolkit library is                                                     *
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
* $Id$                  *
*********************************************************************************/

/*
  Notes:
    Please report errors to: develop@nvutiec.com
 
    The widget sends a SEL_CHANGED width:
    handle(this,FXSEL(SEL_CHANGED,selector),&actualDate);
    to its parent if onCmdPrev,onCmdNext,onCmdToday is activated
    If a date is selected width the mouse, it sends 
    handle(this,FXSEL(SEL_LEFTBUTTONPRESS,selector),&selectedDate);
    to its parent;
 
 
 
  ToDo:
  
    * extra FXArrowButtons for month and year;
    * nicer "today" circle :)
    * colormanagement
*/




#ifndef FXDATATIMEPICKER_H
#define FXDATATIMEPICKER_H

#ifndef FXFRAME_H
#include "fx.h"
#endif
#include "FXDate.h"
namespace FX
  {


  class FXDatePickerLabel : public FXLabel
    {
      FXDECLARE(FXDatePickerLabel)
    private:
      FXint selector;
      FXObject* target;
      FXbool istoday;
      FXColor todayColor;
      FXint dayvalue;
    protected:
      FXDatePickerLabel()
      {
        selector=0;
        target=NULL;
        istoday=false;
        dayvalue=-1;
      };
    private:
      FXDatePickerLabel(const FXDatePickerLabel&);
      FXDatePickerLabel &operator=(const FXDatePickerLabel&);

    public:
      /* Message Callbacks */
      long onClick(FXObject*,FXSelector,void*);
      long onLeftBtnPress(FXObject*,FXSelector,void*);
      long onPaint(FXObject*,FXSelector,void*);
    public:
      enum {
        ID_NOP=FXLabel::ID_LAST,ID_LAST
      };
    public:

      /// Construct
      FXDatePickerLabel(FXComposite* p,FXObject* tgt=NULL,FXSelector sel=0,const FXString& text="",FXIcon* ic=0,FXuint opts=LABEL_NORMAL,FXint x=0,FXint y=0,FXint w=0,FXint h=0,FXint pl=DEFAULT_PAD,FXint pr=DEFAULT_PAD,FXint pt=DEFAULT_PAD,FXint pb=DEFAULT_PAD);
      // This Label is today
      void setToday(FXbool value = true)
      {
        istoday=value;
      }
      FXbool isToday()
      {
        return istoday;
      }

      // Setting color for today
      void setTodayColor(FXColor c)
      {
        todayColor = c;
      }
      // Getting color for today
      FXColor getTodayColor()
      {
        return todayColor;
      }
      // setting value for day
      void setDayValue(FXint dv)
      {
        dayvalue=dv;
      }
      // getting value for day
      FXint getDayValue()
      {
        return dayvalue;
      }
    };

  /*
   
    FXDatePickerHoliday has one method
    
    bool isHoliday(const FXDate &date);
    return true .. the given day is a holliday
           false .. the given day no a holliday
   
    for coloring holidays :)
   
   */
  class FXDatePickerHoliday
    {
    public:
      FXDatePickerHoliday()
      {}
      ;

      virtual FXbool isHoliday(FXDate &date) = 0;

    };

  class FXDatePicker : public FXVerticalFrame
    {
      FXDECLARE(FXDatePicker)
    private:

      FXint selector;
      FXObject* target;
      FXLabel *monthname, *year;
      FXColor monthbkcolor,yearbkcolor;
      FXColor selectColor,offdaycolor, holidaycolor;
      FXButton * next,*prev;
      FXButton *nextYear, *prevYear;
      FXHorizontalFrame * navframe;
      FXObjectList dayfields;
      FXColor boderColorToday;
      FXDatePickerLabel * selected;
      FXDate actualDate;
      FXDate selectedDate; // current selected day
      FXDatePickerHoliday *holiday;    /* will be deleted automaticly */


    protected:
      FXDatePicker();
    private:
      FXDatePicker(const FXDatePicker&);
      FXDatePicker &operator=(const FXDatePicker&);
      ~FXDatePicker();

    public:
      /* Message Callbacks */
      long onLabelSelected(FXObject*,FXSelector,void*);
      long onCmdPrevYear(FXObject*,FXSelector,void*);
      long onCmdPrev(FXObject*,FXSelector,void*);
      long onCmdNext(FXObject*,FXSelector,void*);
      long onCmdNextYear(FXObject*,FXSelector,void*);
      long onCmdToday(FXObject*,FXSelector,void*);


    public:
      enum {
        ID_NOP=FXVerticalFrame::ID_LAST,ID_BUTTON_PREV,ID_BUTTON_NEXT,ID_BUTTON_TODAY,ID_LABEL,
        ID_BUTTON_PREVY, ID_BUTTON_NEXTY, ID_LAST
      };
    public:

      /// Construct color well with initial color clr
      FXDatePicker(FXComposite* p,FXDatePickerHoliday *holi=NULL,FXObject* tgt=NULL,FXSelector sel=0,FXuint opts=0,FXint x=0,FXint y=0,FXint w=0,FXint h=0,FXint pl=DEFAULT_PAD,FXint pr=DEFAULT_PAD,FXint pt=DEFAULT_PAD,FXint pb=DEFAULT_PAD);

      /// Create server-side resources
      virtual void create();

      /// Detach server-side resources
      virtual void detach();

      void fill(FXDate &adate);
      // reset picker to today
      void reset();
      // getting selected FXDate
      FXDate getSelectedDate();
      // setting selected FXDate   //
      FXDate setSelectedDate(FXDate d); //

      // setting country dependend object to check holidays
      void setHolidayObject(FXDatePickerHoliday *obj)
      {
        holiday=obj;
      }

      // setting colors
      void setMonthBkColor(FXColor col);
      void setYearBkColor(FXColor col);
      void setSelectionColor(FXColor col);
      void setOffdayColor(FXColor col);
      void setHolidayColor(FXColor col);

      // getting colors
      FXColor getMonthBkColor()
      {
        return  monthbkcolor;
      }
      FXColor getYearBkColor()
      {
        return  yearbkcolor;
      }
      FXColor getSelectionColor()
      {
        return  selectColor;
      }
      FXColor getOffdayColor()
      {
        return  offdaycolor;
      }
      FXColor getHolidayColor()
      {
        return  holidaycolor;
      }

    };

}
#endif

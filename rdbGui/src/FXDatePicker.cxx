/********************************************************************************
 * Author:  Christoph Singewald                                                 *
 * develop@chenos.com                                                           *
 * www.chenos.com                                                               *
 * see FXDatePicker.h                                                           *
 ********************************************************************************/



#include "FXDatePicker.h"
#include "Icons.h"

using namespace FX;

/*******************************************************************************/

namespace FX
  {


  FXDEFMAP(FXDatePickerLabel) FXDatePickerLabelMap[]={
        FXMAPFUNC(SEL_CLICKED,0,FXDatePickerLabel::onClick),
        FXMAPFUNC(SEL_LEFTBUTTONPRESS,0,FXDatePickerLabel::onLeftBtnPress),
        FXMAPFUNC(SEL_PAINT,0,FXDatePickerLabel::onPaint)

      };

  // Object implementation
  FXIMPLEMENT(FXDatePickerLabel,FXLabel,FXDatePickerLabelMap,ARRAYNUMBER(FXDatePickerLabelMap))

  FXDatePickerLabel::FXDatePickerLabel(FXComposite* p,FXObject* tgt,FXSelector sel,const FXString& text,FXIcon* ic,FXuint opts,FXint x,FXint y,FXint w,FXint h,FXint pl,FXint pr,FXint pt,FXint pb) :
      FXLabel(p,text,ic,opts,x,y,w,h,pl,pr,pt,pb)
  {
    selector = sel;
    target = tgt;
    istoday = false;
    dayvalue=-1;

    /** DefaultColors **/
    todayColor =FXRGB(200,0,0);
  }

  long FXDatePickerLabel::onClick(FXObject*, FXSelector, void*)
  {

    return 0;
  }

  long FXDatePickerLabel::onLeftBtnPress(FXObject*,FXSelector,void*)
  {

    if(target && target->handle(this,FXSEL(SEL_LEFTBUTTONPRESS,selector),0))
      {
        return 1;
      }
    return 0;
  }

  long FXDatePickerLabel::onPaint(FXObject* obj,FXSelector sel,void* ptr)
  {

    FXLabel::onPaint(obj,sel,ptr);
    FXEvent   *ev=(FXEvent*)ptr;
    FXDCWindow dc(this,ev);
    if(istoday)
      {
        dc.setForeground(todayColor);
        dc.drawRectangle (0,0, getWidth()-1, getHeight()-1);
      }

    return 1;
  }

  /*******************************************************************************/


  FXDEFMAP(FXDatePicker) FXDatePickerMap[]={
        FXMAPFUNC(SEL_LEFTBUTTONPRESS,FXDatePicker::ID_LABEL,FXDatePicker::onLabelSelected),
        FXMAPFUNC(SEL_COMMAND,FXDatePicker::ID_BUTTON_PREVY,FXDatePicker::onCmdPrevYear),
        FXMAPFUNC(SEL_COMMAND,FXDatePicker::ID_BUTTON_PREV,FXDatePicker::onCmdPrev),
        FXMAPFUNC(SEL_COMMAND,FXDatePicker::ID_BUTTON_TODAY,FXDatePicker::onCmdToday),
        FXMAPFUNC(SEL_COMMAND,FXDatePicker::ID_BUTTON_NEXT,FXDatePicker::onCmdNext),
        FXMAPFUNC(SEL_COMMAND,FXDatePicker::ID_BUTTON_NEXTY,FXDatePicker::onCmdNextYear)
      };

  // Object implementation
  FXIMPLEMENT(FXDatePicker,FXVerticalFrame,FXDatePickerMap,ARRAYNUMBER(FXDatePickerMap))



  /*******************************************************************************/


  // Init
  FXDatePicker::FXDatePicker()
  {
    selected=NULL;
  }



  FXDatePicker::FXDatePicker(FXComposite* p,FXDatePickerHoliday *holi,FXObject* tgt,FXSelector sel,FXuint opts,FXint x,FXint y,FXint w,FXint h,FXint pl,FXint pr,FXint pt,FXint pb):
      FXVerticalFrame(p,opts,x,y,w,h,pl,pr,pt,pb)
  {

    holiday = holi;
    selector = sel;
    target =tgt;
    selected =NULL;
    /* Default Colors */
    monthbkcolor = FXRGB(184,186,222);
    yearbkcolor = FXRGB(184,186,222);
    boderColorToday = FXRGB(200,0,0);
    selectColor = FXRGB(200,200,200);
    offdaycolor = FXRGB(255,0,0);
    holidaycolor = FXRGB(0,0,255);


    FXFont * f = getApp()->getNormalFont();
    FXint textheight = f->getTextHeight ("><");
    FXint textwidth = f->getTextWidth("30"); // possible max width
    
    FXGIFIcon *dArrowLeftIc = new FXGIFIcon(getApp(), darrowLeft);
    FXGIFIcon *arrowLeftIc = new FXGIFIcon(getApp(), arrowLeft);
    FXGIFIcon *arrowRightIc = new FXGIFIcon(getApp(), arrowRight);
    FXGIFIcon *dArrowRightIc = new FXGIFIcon(getApp(), darrowRight);

    setBackColor(FXRGB(255,255,255));
    /** Creating Components **/
    navframe = new FXHorizontalFrame(this,FRAME_GROOVE|LAYOUT_FILL_X,0,0,textheight,0,0,0,0,0);
    navframe->setBackColor(yearbkcolor);
    prevYear = new FXButton(navframe,"",dArrowLeftIc, this, ID_BUTTON_PREVY,FRAME_THICK|FRAME_RAISED|LAYOUT_CENTER_Y,0,0,0,0,0,0,0,0);
    prev = new FXButton(navframe,"",arrowLeftIc,this, ID_BUTTON_PREV,FRAME_THICK|FRAME_RAISED|LAYOUT_CENTER_Y,0,0,0,0,0,0,0,0);
    monthname = new FXLabel(navframe,"Januar",NULL,FRAME_NONE|LAYOUT_FILL_X);
    monthname->setBackColor(monthbkcolor);
    year = new FXLabel(navframe,"2004",NULL,FRAME_NONE);
    year->setBackColor(yearbkcolor);
    next = new FXButton(navframe,"",arrowRightIc, this, ID_BUTTON_NEXT,FRAME_THICK|FRAME_RAISED| LAYOUT_CENTER_Y,0,0,0,0,0,0,0,0);
    nextYear = new FXButton(navframe,"", dArrowRightIc,this, ID_BUTTON_NEXTY, FRAME_THICK|FRAME_RAISED| LAYOUT_CENTER_Y| LAYOUT_RIGHT,0,0,0,0,0,0,0,0);
    FXMatrix * days =new FXMatrix(this,7,MATRIX_BY_COLUMNS| LAYOUT_FILL_COLUMN| LAYOUT_SIDE_TOP| LAYOUT_SIDE_RIGHT|LAYOUT_FILL_X,0,0,0,0,0,0,0,0);
    days->setBackColor(FXRGB(255,255,255));


    FXLabel * l;
    l = new FXLabel(days,FXDate::getShortDayName(1),NULL,FRAME_NONE);
    l->setBackColor(FXRGB(255,255,255));
    l = new FXLabel(days,FXDate::getShortDayName(2),NULL,FRAME_NONE);
    l->setBackColor(FXRGB(255,255,255));
    l = new FXLabel(days,FXDate::getShortDayName(3),NULL,FRAME_NONE);
    l->setBackColor(FXRGB(255,255,255));
    l = new FXLabel(days,FXDate::getShortDayName(4),NULL,FRAME_NONE);
    l->setBackColor(FXRGB(255,255,255));
    l = new FXLabel(days,FXDate::getShortDayName(5),NULL,FRAME_NONE);
    l->setBackColor(FXRGB(255,255,255));
    l = new FXLabel(days,FXDate::getShortDayName(6),NULL,FRAME_NONE);
    l->setBackColor(FXRGB(255,255,255));
    l->setTextColor(offdaycolor);
    l = new FXLabel(days,FXDate::getShortDayName(0),NULL,FRAME_NONE);
    l->setBackColor(FXRGB(255,255,255));
    l->setTextColor(offdaycolor);

    FXint count  = 0;
    for(FXint i=0; i<42;i++)
      {
        FXDatePickerLabel *l = new FXDatePickerLabel(days,this,ID_LABEL," ",NULL,FRAME_NONE | JUSTIFY_CENTER_X | JUSTIFY_CENTER_Y | LAYOUT_FILL_X| LAYOUT_FILL_Y |LAYOUT_CENTER_X | LAYOUT_CENTER_Y,0,0,textwidth,0);
        l->setBackColor(FXRGB(255,255,255));
        dayfields.append(l);
        if(count == 5 || count==6)
          {
            l->setTextColor(offdaycolor);
          }
        if(count==6)
          count = 0;
        else
          count++;

      }

    selectedDate=actualDate;
    fill(actualDate);
    prevYear->setWidth(prev->getWidth());
    prevYear->setHeight(prev->getHeight());
  }

  FXDatePicker::~FXDatePicker()
  {
    if(holiday)
      delete holiday;
  }



  // Create window
  void FXDatePicker::create()
  {
    FXVerticalFrame::create();
  }


  // Detach window
  void FXDatePicker::detach()
  {
    FXVerticalFrame::detach();

  }
  void FXDatePicker::setMonthBkColor(FXColor col)
  {
    monthbkcolor=col;
    navframe->setBackColor(monthbkcolor);
    monthname->setBackColor(monthbkcolor);
  }
  void FXDatePicker::setYearBkColor(FXColor col)
  {
    yearbkcolor=col;
    navframe->setBackColor(yearbkcolor);
    year->setBackColor(yearbkcolor);
  }
  void FXDatePicker::setSelectionColor(FXColor col)
  {
    selectColor=col;
  }
  void FXDatePicker::setOffdayColor(FXColor col)
  {
    offdaycolor=col;
    FXint count = 0;
    for(FXint i=0; i<42;i++)
      {
        if(count == 5 || count==6)
          {
            ((FXLabel * )dayfields[i])->setTextColor(offdaycolor);
          }
        if(count==6)
          count = 0;
        else
          count++;

      }
  }
  void FXDatePicker::setHolidayColor(FXColor col)
  {
    holidaycolor=col;
  }


  long FXDatePicker::onLabelSelected(FXObject* obj, FXSelector, void*)
  {

    // selected öabel is not valid
    if(((FXDatePickerLabel*)obj)->getDayValue()==-1)
      return 1;

    if(selected)
      {
        selected->setBackColor(getBackColor());
      }
    selected = (FXDatePickerLabel *)obj;
    selected->setBackColor(selectColor);
    // We have a date selected send a message to target
    //FXDate d = getSelectedDate();
    selectedDate.setDay(selected->getDayValue());
    selectedDate.setMonth(actualDate.getMonth());
    selectedDate.setYear(actualDate.getYear());
    if(target)
      target->handle(this,FXSEL(SEL_LEFTBUTTONPRESS,selector),&selectedDate);

    return 1;

    //return 0;
  }
  
  long FXDatePicker::onCmdPrevYear(FXObject*, FXSelector, void*)
  {
    actualDate.setYear(actualDate.getYear()-1);
    fill(actualDate);
    if(target)
      target->handle(this,FXSEL(SEL_CHANGED,selector),&actualDate);
    return 0;
  }
  
  long FXDatePicker::onCmdPrev(FXObject*, FXSelector, void*)
  {
    actualDate.setMonth(actualDate.getMonth()-1);
    fill(actualDate);
    if(target)
      target->handle(this,FXSEL(SEL_CHANGED,selector),&actualDate);
    return 0;
  }
  long FXDatePicker::onCmdToday(FXObject*, FXSelector, void*)
  {

    actualDate = FXDate();
    fill(actualDate);
    if(target)
      target->handle(this,FXSEL(SEL_CHANGED,selector),&actualDate);

    return 0;
  }
  long FXDatePicker::onCmdNext(FXObject*, FXSelector, void*)
  {
    actualDate.setMonth(actualDate.getMonth()+1);
    fill(actualDate);
    if(target)
      target->handle(this,FXSEL(SEL_CHANGED,selector),&actualDate);
    return 0;
  }
  
  long FXDatePicker::onCmdNextYear(FXObject*, FXSelector, void*)
  {
    actualDate.setYear(actualDate.getYear()+1);
    fill(actualDate);
    if(target)
      target->handle(this,FXSEL(SEL_CHANGED,selector),&actualDate);
    return 0;
  }

  
  void FXDatePicker::reset()
  {
    actualDate = FXDate();
    fill(actualDate);

  }
  void FXDatePicker::fill(FXDate &adate)
  {

    FXint i = 0;
    FXString buffer ="";
    FXDate firstdayofmonth(adate.getYear(),adate.getMonth(),1);
    FXDate today;
    FXDate tempdate =firstdayofmonth;
    if (selected)
      {
        selected->setBackColor(getBackColor());
      }
    selected=NULL;

    monthname->setText(adate.getMonthname());
    buffer.format("%04d",adate.getYear());
    year->setText(buffer);

    FXint dayindex = (firstdayofmonth.getDayOfWeek()==0)?6:firstdayofmonth.getDayOfWeek()-1;

    // Delete unused

    for(i=0; i<dayfields.no();i++)
      {
        ((FXLabel * )dayfields[i])->setText(" ");
        ((FXDatePickerLabel * )dayfields[i])->setToday(false);
        ((FXDatePickerLabel * )dayfields[i])->setDayValue(-1);
      }


    FXDate d=adate;
    for(i=dayindex; i<adate.getDaysOfMonth(adate.getMonth())+dayindex;i++)
      {
        buffer.format("%d",i-dayindex+1);
        ((FXLabel * )dayfields[i])->setText(buffer);
        ((FXDatePickerLabel * )dayfields[i])->setDayValue(i-dayindex+1);

        d.setDay(i-dayindex+1);
        if (selectedDate == d)
          {
            selected = (FXDatePickerLabel * )dayfields[i];
            selected->setBackColor(selectColor);
          }
        if(today == d)
          {
            ((FXDatePickerLabel * )dayfields[i])->setToday();
            ((FXDatePickerLabel * )dayfields[i])->setTodayColor(boderColorToday);
          }
        if(holiday)
          {
            tempdate.setDay(i-dayindex+1);
            if(holiday->isHoliday(tempdate))
              {
                ((FXLabel * )dayfields[i])->setTextColor(holidaycolor);
              }
            else if(tempdate.getDayOfWeek()==0 || tempdate.getDayOfWeek()==6)
              {
                ((FXLabel * )dayfields[i])->setTextColor(offdaycolor);
              }
            else
              ((FXLabel * )dayfields[i])->setTextColor(FXRGB(0,0,0));
          }



      }
  }

  FXDate  FXDatePicker::getSelectedDate()
  {
    /*if(selected)
      {
      actualDate.setDay(selected->getDayValue());
      }
      return actualDate;*/
    return selectedDate;

  }

  FXDate  FXDatePicker::setSelectedDate(FXDate d)
  {
    actualDate=d;
    selectedDate=d;
    fill(actualDate);
    return actualDate;
  }


}

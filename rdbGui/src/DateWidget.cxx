#include "FXDatePicker.h"
#include "DateWidget.h"

/// This dialog is used to show a nice calendar
/// to choose a date and time
class CalDialog : public FXDialogBox
{
  FXDECLARE(CalDialog)

// Services
public:
  CalDialog(FXApp* app);
  ~CalDialog() { }

  long onCmdAccept(FXObject* sender,FXSelector sel,void* ptr);

  FXDate getDate(){return m_date;};

private:
  FXDatePicker * dp;  
  FXDate m_date;
  CalDialog() {}   
};

FXDEFMAP(CalDialog) CalDialogMap[]={
    FXMAPFUNC(SEL_COMMAND, FXDialogBox::ID_ACCEPT, CalDialog::onCmdAccept)    
};

FXIMPLEMENT(CalDialog, FXDialogBox, CalDialogMap,ARRAYNUMBER(CalDialogMap))


CalDialog::CalDialog(FXApp* app)
: FXDialogBox(app, "Choose a date and time", DECOR_TITLE|DECOR_BORDER)
{
  // Mainframe for deco
  FXVerticalFrame *mainframe =new FXVerticalFrame(this, FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y|LAYOUT_TOP|LAYOUT_LEFT);


  dp = new FXDatePicker(mainframe);
//  dp->setMonthBkColor(FXRGB(158,234,231));
//  dp->setYearBkColor(FXRGB(158,234,231));

  dp->setMonthBkColor(FXRGB(255,204,49));
  dp->setYearBkColor(FXRGB(255,204,49));
  
  FXHorizontalFrame *buttonframe =new FXHorizontalFrame(mainframe,LAYOUT_FILL_Y| LAYOUT_FILL_X|LAYOUT_TOP|LAYOUT_RIGHT,0,0,0,0,0,0,0,0);
  new FXButton(buttonframe ,"&OK",NULL,this,FXDialogBox::ID_ACCEPT,FRAME_RAISED|FRAME_THICK|LAYOUT_RIGHT|LAYOUT_CENTER_Y|LAYOUT_FIX_WIDTH|LAYOUT_FIX_HEIGHT,0,0,70,25);
 
}


long  CalDialog::onCmdAccept(FXObject* sender,FXSelector sel,void* ptr)
{
    FXDialogBox::onCmdAccept(sender,sel,ptr);
    m_date = dp->getSelectedDate();

    return 1;
}



class DateField : public FXHorizontalFrame
{
  FXDECLARE(DateField)

// Services
public:
   enum{
    ID_LIST=FXFrame::ID_LAST,
    ID_SET
    };

  DateField(FXComposite* owner);
  ~DateField() { }

  long onCmdSet(FXObject* sender,FXSelector sel,void* ptr);
  void setDate(std::string d){m_dateField->setText(d.c_str());};
  std::string getDate(){return m_dateField->getText().text();};

private:
  FXTextField* m_dateField;
  FXButton* m_startCal;
  DateField() {}   
};

FXDEFMAP(DateField) DateFieldMap[]={
    FXMAPFUNC(SEL_COMMAND, DateField::ID_SET, DateField::onCmdSet)    
};

FXIMPLEMENT(DateField, FXHorizontalFrame, DateFieldMap,ARRAYNUMBER(DateFieldMap))

DateField::DateField(FXComposite* owner)
: FXHorizontalFrame(owner, LAYOUT_FILL_X, LAYOUT_FILL_Y, 0, 0, 0, 0, 0, 0, 0, 0)
{
  m_dateField = new FXTextField((FXComposite*)this,3,NULL,0,TEXTFIELD_NORMAL|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
  m_startCal = new FXButton((FXComposite*)this, "...", NULL, this, ID_SET);
  
  m_dateField->setText("");
}


long  DateField::onCmdSet(FXObject* sender,FXSelector sel,void* ptr)
{

  CalDialog* test = new CalDialog(getApp());
  test->show();
 
  if (test->execute(PLACEMENT_OWNER) != 0)
  {
    FXDate temp = test->getDate();

    char name[50];
    sprintf(name,"%d-%d-%d %.2d:%.2d:%.2d", 
    temp.getYear(),
    temp.getMonth(),
    temp.getDay(),
    temp.getHour(),
    temp.getMinute(),
    temp.getSecond());
    std::string res(name);
    m_dateField->setText(name);
    
    return 1; 
  }

  return 1;
}

DateWidget::DateWidget(FXComposite* parent, rdbModel::Column *column)
{
  m_column = column;
  m_widget = new FXHorizontalFrame(parent,LAYOUT_FILL_X|LAYOUT_FILL_COLUMN,
      0, 0, 0, 0, 0, 0, 0, 0);
 
  m_dateField = new DateField(m_widget);
  m_dateField->setVSpacing(0);
  m_dateField->setHSpacing(0);
  m_dateField->setBackColor(parent->getBackColor());
      
  m_widget->create();
}

std::string DateWidget::getValue()
{
  return m_dateField->getDate();
}

void DateWidget::setValue(std::string val)
{
    m_dateField->setDate(val);
}  

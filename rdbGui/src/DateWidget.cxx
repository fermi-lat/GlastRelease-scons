#include "DateWidget.h"

DateWidget::DateWidget(FXComposite* parent, rdbModel::Column *column)
{
  m_column = column;
  m_widget = new FXHorizontalFrame(parent);

  m_year = new FXComboBox(m_widget,4,NULL,0,
			  COMBOBOX_STATIC|FRAME_SUNKEN|FRAME_THICK|LAYOUT_SIDE_TOP);
  m_year->setNumVisible(10);
  for(int i=0; i<68; i++){
    int j = 1970;
    char name[50];
    sprintf(name,"%04d",j+i);
    m_year->appendItem(name);
  }

  m_month = new FXComboBox(m_widget,2,NULL,0,
			   COMBOBOX_STATIC|FRAME_SUNKEN|FRAME_THICK|LAYOUT_SIDE_TOP);
  m_month->setNumVisible(12);
  for(int i=0; i<12; i++){
    char name[50];
    sprintf(name,"%02d",i+1);
    m_month->appendItem(name);
  }

  m_day = new FXComboBox(m_widget,2,NULL,0,
			 COMBOBOX_STATIC|FRAME_SUNKEN|FRAME_THICK|LAYOUT_SIDE_TOP);
  m_day->setNumVisible(10);
  for(int i=0; i<31; i++){
    char name[50];
    sprintf(name,"%02d",i+1);
    m_day->appendItem(name);
  }
}

std::string DateWidget::getValue()
{
  char name[50];
  sprintf(name,"%s-%s-%s", 
	  m_year->getText().text(),
	  m_month->getText().text(),
	  m_day->getText().text());
  std::string res(name);
  return res;
}

void DateWidget::setValue(std::string val)
{
  int year, month, day;

  sscanf(val.c_str(),"%d-%d-%d", &year, &month, &day);  

  m_year->setText(FXStringVal(year));
  m_month->setText(FXStringVal(month));
  m_day->setText(FXStringVal(day));
}  

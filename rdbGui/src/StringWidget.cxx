#include "StringWidget.h"

StringWidget::StringWidget(FXComposite* parent, rdbModel::Column *column)
{
  m_column = column;
  m_widget = new FXTextField((FXComposite*)parent,3,NULL,0,TEXTFIELD_NORMAL|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
}

std::string StringWidget::getValue()
{
  std::string res((m_widget->getText()).text());
  return res;
}

void StringWidget::setValue(std::string val)
{
  m_widget->setText(val.c_str());
}  

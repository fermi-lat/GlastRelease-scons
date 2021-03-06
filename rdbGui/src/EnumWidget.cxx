#include "EnumWidget.h"
#include "rdbModel/Tables/Datatype.h"

EnumWidget::EnumWidget(FXComposite* parent, rdbModel::Column *column)
{
  rdbModel::Datatype* dt = column->getDatatype();
  rdbModel::Enum* en = dt->getEnum();

  m_column = column;
  m_widget = new FXComboBox(parent,0, NULL, 0,FRAME_SUNKEN|FRAME_THICK|LAYOUT_SIDE_TOP|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
  m_widget->setNumVisible(en->getChoices().size());
  if (en->choicesRequired())
    m_widget->setEditable(FALSE); 	
  m_widget->create();
  
  for(unsigned int i=0; i<en->getChoices().size(); i++){
    m_widget->appendItem((en->getChoices()[i]).c_str());
  }
}

std::string EnumWidget::getValue()
{
  std::string res(m_widget->getText().text());
  return res;
}

void EnumWidget::setValue(std::string val)
{ 
  m_widget->setText(FXString(val.c_str()));
}  

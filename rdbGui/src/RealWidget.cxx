#include "RealWidget.h"
#include "rdbModel/Tables/Datatype.h"

RealWidget::RealWidget(FXComposite* parent, rdbModel::Column *column)
{
  rdbModel::Datatype* dt = column->getDatatype();

  m_column = column;

  m_widget = new FXRealSpinner(parent,20,
                               NULL,0,
                               SPIN_NORMAL|FRAME_SUNKEN|FRAME_THICK|LAYOUT_SIDE_TOP);
//spinner->setRange(1,20);

}

std::string RealWidget::getValue()
{
  char name[50];
  sprintf(name,"%f", m_widget->getValue());
  std::string res(name);
  return res;
}

void RealWidget::setValue(std::string val)
{
  float res;
  sscanf(val.c_str(),"%f", &res);
  m_widget->setValue(res);
}  

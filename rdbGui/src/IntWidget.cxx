#include "IntWidget.h"
#include "rdbModel/Tables/Datatype.h"
#include <limits.h>

IntWidget::IntWidget(FXComposite* parent, rdbModel::Column *column)
{
  rdbModel::Datatype* dt = column->getDatatype();

  m_column = column;

  m_widget = new FXSpinner(parent,0,
                           NULL,0,
                           SPIN_NORMAL|FRAME_SUNKEN|FRAME_THICK|LAYOUT_SIDE_TOP|
                           LAYOUT_FILL_X);
  m_widget->create();


  /// Lets now set some possible numeric restrictions
  switch(dt->getRestrict()){  
    case rdbModel::Datatype::RESTRICTnonneg :
      {
        //XXX to change with a more proper limit
        m_widget->setRange(0,INT_MAX);
        break;
      };
    case rdbModel::Datatype::RESTRICTpos :
      {
        //XXX to change with a more proper limit
        m_widget->setRange(1,INT_MAX);
        break;
      };
    case rdbModel::Datatype::RESTRICTinterval :
      {
        //XXX to change with a more proper limit
        std::string min, max;
        int mi, ma;
        dt->getInterval(min, max);
        sscanf(min.c_str(),"%d", &mi);     
        sscanf(max.c_str(),"%d", &ma);     
        m_widget->setRange(mi,ma);
        break;
      };
  }

}

std::string IntWidget::getValue()
{
  char name[50];
  sprintf(name,"%d", m_widget->getValue());
  std::string res(name);
  return res;
}

void IntWidget::setValue(std::string val)
{
  int res;
  sscanf(val.c_str(),"%d", &res);
  m_widget->setValue(res);
}  

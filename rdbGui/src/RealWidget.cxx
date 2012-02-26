#include <limits.h>
#include "RealWidget.h"
#include "rdbModel/Tables/Datatype.h"

RealWidget::RealWidget(FXComposite* parent, rdbModel::Column *column)
{
  rdbModel::Datatype* dt = column->getDatatype();

  m_column = column;

  m_widget = new FXRealSpinner(parent,0,
                               NULL,0,
                               SPIN_NORMAL|FRAME_SUNKEN|FRAME_THICK|LAYOUT_SIDE_TOP|
                               LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
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
    default:
      {};
  }

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

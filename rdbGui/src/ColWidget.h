#ifndef COLWIDGET_H
#define COLWIDGET_H

#include <iostream>
#include "fx.h"
#include "rdbModel/Tables/Column.h"

class ColWidget{
 public:
  rdbModel::Column* getColumn(){return m_column;};
  virtual std::string getValue() =  0;
  virtual void setValue(std::string) = 0;
  virtual FXWindow* getWidget() = 0;

 protected:
  rdbModel::Column* m_column;
};

#endif

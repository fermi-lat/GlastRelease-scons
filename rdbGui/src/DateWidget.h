#ifndef DATEWIDGET_H
#define DATEWIDGET_H

#include "ColWidget.h"

class DateField;

class DateWidget : public ColWidget{
  public:
  DateWidget(FXComposite*, rdbModel::Column*);
  virtual std::string getValue();
  virtual void setValue(std::string);  
  virtual FXWindow* getWidget(){return m_widget;} 
  
 private:
  FXHorizontalFrame* m_widget;  
  DateField* m_dateField;
};

#endif

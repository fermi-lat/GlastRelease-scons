#ifndef STRINGWIDGET_H
#define STRINGWIDGET_H
#include "ColWidget.h"

class StringWidget : public ColWidget{
 public:
  StringWidget(FXComposite*, rdbModel::Column*);
  virtual std::string getValue();
  virtual void setValue(std::string);  
  virtual FXWindow* getWidget(){return m_widget;} 
  
 private:
  FXTextField* m_widget;  
};

#endif

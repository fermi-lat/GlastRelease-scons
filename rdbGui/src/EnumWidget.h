#ifndef ENUMWIDGET_H
#define ENUMWIDGET_H
#include "ColWidget.h"

class EnumWidget : public ColWidget{
 public:
  EnumWidget(FXComposite*, rdbModel::Column*);
  virtual std::string getValue();
  virtual void setValue(std::string);  
  virtual FXWindow* getWidget(){return m_widget;} 
  
 private:
  FXComboBox* m_widget;  
};


#endif

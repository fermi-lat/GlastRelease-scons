#ifndef REALWIDGET_H
#define REALWIDGET_H

#include "ColWidget.h"

class RealWidget : public ColWidget{
 public:
  RealWidget(FXComposite*, rdbModel::Column*);
  virtual std::string getValue();
  virtual void setValue(std::string);  
  virtual FXWindow* getWidget(){return m_widget;} 
  
 private:
  FXRealSpinner* m_widget;  
};



#endif

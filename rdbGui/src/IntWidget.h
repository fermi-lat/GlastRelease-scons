#ifndef INTWIDGET_H
#define INTWIDGET_H

#include "ColWidget.h"

class IntWidget : public ColWidget{
 public:
  IntWidget(FXComposite*, rdbModel::Column*);
  virtual std::string getValue();
  virtual void setValue(std::string);  
  virtual FXWindow* getWidget(){return m_widget;} 
  
 private:
  FXSpinner* m_widget;  
};


#endif

#ifndef COLWIDGETFACTORY_H
#define COLWIDGETFACTORY_H

#include "ColWidget.h"

class ColWidgetFactory{
 public:
  ColWidget* createStringWidget(FXComposite*, rdbModel::Column*);
  ColWidget* createDateWidget(FXComposite*, rdbModel::Column*);
  ColWidget* createEnumWidget(FXComposite*, rdbModel::Column*);
  ColWidget* createRealWidget(FXComposite*, rdbModel::Column*);
  ColWidget* createIntWidget(FXComposite*, rdbModel::Column*);
  ColWidget* createColWidget(FXComposite*, rdbModel::Column*);
};

#endif

#include "ColWidgetFactory.h"
#include "StringWidget.h"
#include "DateWidget.h"
#include "EnumWidget.h"
#include "RealWidget.h"
#include "IntWidget.h"

ColWidget* ColWidgetFactory::createStringWidget(FXComposite* parent,rdbModel::Column *column)
{
  return new StringWidget(parent, column);
}

ColWidget* ColWidgetFactory::createDateWidget(FXComposite* parent, rdbModel::Column *column)
{
  return new DateWidget(parent, column);
}

ColWidget* ColWidgetFactory::createEnumWidget(FXComposite* parent, rdbModel::Column *column)
{
  return new EnumWidget(parent, column);
}

ColWidget* ColWidgetFactory::createRealWidget(FXComposite* parent, rdbModel::Column *column)
{
  return new RealWidget(parent, column);
}

ColWidget* ColWidgetFactory::createIntWidget(FXComposite* parent, rdbModel::Column *column)
{
  return new IntWidget(parent, column);
}

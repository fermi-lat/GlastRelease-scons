#ifndef INSERTDIALOG_H
#define INSERTDIALOG_H

#include "fx.h"
#include <vector>
#include "rdbModel/Management/Visitor.h"
#include "rdbModel/Rdb.h"
#include "rdbModel/Tables/Table.h"
#include "rdbModel/Tables/Assertion.h"
#include "rdbModel/Tables/Column.h"

class ColWidgetFactory;
class ColWidget;

class InsertDialog: public FXDialogBox,public rdbModel::Visitor
{
  FXDECLARE(InsertDialog)
 public:
 
   enum{
    ID_LIST=FXDialogBox::ID_LAST,
    ID_GO
   };

 
  InsertDialog(FXWindow *owner);
  
  long onGoPress(FXObject*,FXSelector,void*);
  
 protected:
  InsertDialog(){}
  InsertDialog(const InsertDialog&){}
  
 private:

  // This is the main vertical frame of the insert form
  FXVerticalFrame *m_uiRpanel;

  // This is the widget matrix with the insert fields
  FXMatrix *m_matrix;
  ColWidgetFactory* m_factory;  

  /// Visitors to build the insert widgets  
  rdbModel::Visitor::VisitorState visitRdb(rdbModel::Rdb *rdb);
  rdbModel::Visitor::VisitorState visitTable(rdbModel::Table *table);
  rdbModel::Visitor::VisitorState visitColumn(rdbModel::Column *column);
  rdbModel::Visitor::VisitorState visitIndex(rdbModel::Index *index);
  rdbModel::Visitor::VisitorState visitAssertion(rdbModel::Assertion *assertion);

  // The vector of widgets for the insert operation
  std::vector<ColWidget*> m_widgets;
};


#endif

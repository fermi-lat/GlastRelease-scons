#ifndef INSERTDIALOG_H
#define INSERTDIALOG_H

#include "fx.h"
#include <vector>
#include "rdbModel/Management/Visitor.h"
#include "rdbModel/Rdb.h"
#include "rdbModel/Tables/Table.h"
#include "rdbModel/Tables/Assertion.h"
#include "rdbModel/Tables/Column.h"
#include "rdbModel/Db/Connection.h"
#include "rdbModel/Db/ResultHandle.h"

class ColWidgetFactory;
class ColWidget;
class LogText;

class InsertDialog: public FXDialogBox,public rdbModel::Visitor
{
  FXDECLARE(InsertDialog)
 public:
 
   enum{
    ID_LIST=FXDialogBox::ID_LAST,
    ID_GO
   };

 
  InsertDialog(FXApp *owner);
 
  void setConnection(rdbModel::Connection* con){m_connection = con;}
  void setTableName(std::string name){m_tableName = name;} 
  void setUiLog(LogText* log){m_uiLog = log;};
  long onGoPress(FXObject*,FXSelector,void*);

  /// Setter and getter for the last row inserted
  int getLastRow(){return m_lastRow;};
  void setLastRow(int l){m_lastRow = l;};
  
  bool getInsertMode(){return m_insertMode;};
  void setInsertMode(bool m){m_insertMode = m;};

  void fillWithLastRow();

 protected:
  InsertDialog(){}
  InsertDialog(const InsertDialog&){}
  
 private:

  /// This is the main vertical frame of the insert form
  FXVerticalFrame *m_uiRpanel;

  /// This is the widget matrix with the insert fields
  FXMatrix *m_matrix;
  ColWidgetFactory* m_factory;  

  /// Visitors to build the insert widgets  
  rdbModel::Visitor::VisitorState visitRdb(rdbModel::Rdb *rdb);
  rdbModel::Visitor::VisitorState visitTable(rdbModel::Table *table);
  rdbModel::Visitor::VisitorState visitColumn(rdbModel::Column *column);
  rdbModel::Visitor::VisitorState visitIndex(rdbModel::Index *index);
  rdbModel::Visitor::VisitorState visitAssertion(rdbModel::Assertion *assertion);

  
  /// The name of the current table
  std::string m_tableName;
  /// The pointer to the database connection class
  rdbModel::Connection* m_connection;
  /// The vector of widgets for the insert operation
  std::vector<ColWidget*> m_widgets;
  /// The ui log text, to be updated after insertion
  LogText *m_uiLog;
  /// The last row id inserted
  int m_lastRow;                
  /// The mode of this dialog; 1 is insert mode, 0 is update mode
  bool m_insertMode;
  /// The result of the query for the last row inserted 
  rdbModel::ResultHandle* m_result;
};


#endif

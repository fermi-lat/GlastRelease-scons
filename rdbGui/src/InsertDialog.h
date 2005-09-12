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

namespace rdbModel {
  class InsertRow;
  class Supersede;
  class InterRow;
  class Query;
  class Set;
}

class InsertDialog: public FXDialogBox,public rdbModel::Visitor
{
  FXDECLARE(InsertDialog)
 public:
 
   enum{
    ID_LIST=FXDialogBox::ID_LAST,
    ID_GO,
   };

 
  InsertDialog(FXApp *owner);
 
  void setConnection(rdbModel::Connection* con){m_connection = con;}
  void setRdb(rdbModel::Rdb* rdb){m_rdb = rdb;}
  void setTableName(std::string name){m_tableName = name;} 
  void setUiLog(LogText* log){m_uiLog = log;};
  long onGoPress(FXObject*,FXSelector,void*);

  /// Setter and getter for the last table to which a record has been inserted
  std::string getLastTblName(){return m_lastTblName;}
  void setLastTblName(std::string tbName){m_lastTblName=tbName;}
  
  /// Setter and getter for the last row inserted
  int getLastRow(){return m_lastRow;};
  void setLastRow(int l){m_lastRow = l;};
  
  unsigned getInsertMode(){return m_insertMode;};
  void setInsertMode(unsigned m){m_insertMode = m;};

  
  void fillWithRowByKey(std::string ser);
  void fillWithLastRow();

  /// Fills a new row with the last one only in the sticky values ..
  void fillStickyWithLastRow();

 protected:
  InsertDialog(){}
  InsertDialog(const InsertDialog&);
  
 private:

  /// This is the main vertical frame of the insert form
  FXVerticalFrame *m_uiRpanel;

  /// This is the widget matrix with the insert fields
  FXMatrix *m_matrix;
  ColWidgetFactory* m_factory;  

  /// This is the button for multi-insert operation
  FXButton* multi;
    
  /// Visitors to build the insert widgets  
  rdbModel::Visitor::VisitorState visitRdb(rdbModel::Rdb *rdb);
  rdbModel::Visitor::VisitorState visitTable(rdbModel::Table *table);
  rdbModel::Visitor::VisitorState visitColumn(rdbModel::Column *column);
  rdbModel::Visitor::VisitorState visitIndex(rdbModel::Index *index);
  rdbModel::Visitor::VisitorState visitAssertion(rdbModel::Assertion *assertion);
  rdbModel::Visitor::VisitorState visitInsertNew(rdbModel::InsertNew*) {
    return Visitor::VCONTINUE;
  }
  rdbModel::Visitor::VisitorState visitSupersede(rdbModel::Supersede* ) {
    return Visitor::VCONTINUE;
  }
  rdbModel::Visitor::VisitorState visitQuery(rdbModel::Query* ){
    return Visitor::VCONTINUE;
  }
  rdbModel::Visitor::VisitorState visitInterRow(rdbModel::InterRow* ) {
    return Visitor::VCONTINUE;
  }
  rdbModel::Visitor::VisitorState visitSet(rdbModel::Set*) {
    return Visitor::VCONTINUE;
  }

  /// Rdb corresponding to schema (if compatible with connection db)
  rdbModel::Rdb *m_rdb;  
  /// The name of the current table
  std::string m_tableName;
  /// Primary key name (for now compound keys are not considered)
  std::string m_primKey;
  /// The pointer to the database connection class
  rdbModel::Connection* m_connection;
  /// The vector of widgets for the insert operation
  std::vector<ColWidget*> m_widgets;
  /// The ui log text, to be updated after insertion
  LogText *m_uiLog;
  /// Last table increased with a new row
  std::string m_lastTblName;
  /// The last row id inserted
  int m_lastRow;
  /// Key of the selected row
  std::string m_selRow;                
  /// The mode of this dialog; 1 is insert mode, 0 is update mode,
  /// 2 is insertLatest (special kind of insert)
  unsigned m_insertMode;
  /// The result of the query for the last row inserted 
  rdbModel::ResultHandle* m_result;
  /// List of columns set by the service (i.e. by the rdbGui program)
  std::vector<rdbModel::Column* > m_fromService;
};


#endif

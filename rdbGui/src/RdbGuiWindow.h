
#ifndef RDBGUIWINDOW_H
#define RDBGUIWINDOW_H

#include "fx.h"
#include "TableColumnList.h"
#include "QueryFrame.h"
#include "SQLBuffer.h"
#include "LogText.h"
#include "ResultTable.h"
#include "ConnectionDialog.h"

#include "rdbModel/Management/Manager.h"
#include "rdbModel/Management/XercesBuilder.h"
#include "rdbModel/Db/MysqlConnection.h"

#include <cstdio>

class InsertDialog;


// Main Window
class RdbGUIWindow : public FXMainWindow {
  FXDECLARE(RdbGUIWindow)
private:
  FXMenuBar                      *uiMenuBar;               // The menubar
  FXMenuPane                     *uiFilemenu;
  FXMenuPane                     *uiSessmenu;
  FXMenuPane                     *uiActionmenu;
  FXComboBox                     *m_uiDBSelection;         // Showing current connection
  TableColumnList                *uiTblColList;            // List of tables and columns
  QueryFrame                     *searchFrame;             // Frame to perform queries
  SQLBuffer                      *uiEditor;                // SQL editor
  FXbool                          m_shactive;              // Syntax Highlighting set or not
  LogText                        *uiLog;                   // Result log
  ResultTable                    *uiTable;                 // Result Table
  ConnectionDialog               *m_dgNewCon;              // Connection Dialog
  InsertDialog                   *m_dgInsert;              // Insert Dialog
  
  rdbModel::Manager              *m_rdbManager;            // Manager for the rdb
  rdbModel::XercesBuilder        *m_rdbBuilder;            // Builder of the rdb from the xml file
  rdbModel::MysqlConnection      *m_connect;               // Object to connect to a mysql db
  
protected:
  RdbGUIWindow(){}
public:

  // We define additional ID's, starting from the last one used by the base class+1.
  // This way, we know the ID's are all unique for this particular target.
  enum{
    ID_TITLE=FXMainWindow::ID_LAST,
    ID_OVERSTRIKE,
    ID_TREE,
    ID_QUIT,
    ID_SQLEDIT,
    ID_OPENXML,
    ID_OPENCONNECTION,
    ID_CLOSECONNECTION,
    ID_TOGGLERESULT,
    ID_PASTEFROMTABLE,
    ID_TOGGLEAUTOCOMMIT,
    ID_COMMIT,
    ID_ABOUT,
    ID_ROLLBACK, 
    ID_INSERT
    };

  // Message handlers
  long onUpdTitle(FXObject*,FXSelector,void*);
  long onUpdOverstrike(FXObject*,FXSelector,void*);
  long onSQLEdit(FXObject*,FXSelector,void*);
  long onOpenXMLFile(FXObject*,FXSelector,void*);
  long onQuit(FXObject*,FXSelector,void*);
  long onOpenConnection(FXObject*,FXSelector, void*);  
  long onCloseConnection(FXObject*,FXSelector, void*);
  long onQueryFrameUpdate(FXObject *, FXSelector, void*);
  long onInsert(FXObject*,FXSelector, void*);
  // ..eccetera
  
public:

  // RdbGUIWindow constructor
  RdbGUIWindow(FXApp* a);

  // Initialize
  virtual void create();
  
  // RdbGUIWindow destructor
  virtual ~RdbGUIWindow();
  
  
private:

  void loadXMLFile(FXString);
  
  void closeConnection();
  };
  
  
#endif

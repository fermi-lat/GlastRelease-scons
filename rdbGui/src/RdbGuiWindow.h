
#ifndef RDBGUIWINDOW_H
#define RDBGUIWINDOW_H

#include "fx.h"
#include "TableColumnList.h"
#include "QueryFrame.h"
#include "SQLBuffer.h"
#include "LogText.h"
#include "ResultTable.h"
#include "ConnectionDialog.h"

#include "rdbModel/Rdb.h"
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
  FXMenuCommand                  *m_cmdOpenConn;           // Menù item to open a connection
  FXMenuCommand                  *m_cmdCloseConn;          // Menù item to close a connection
  FXMenuPane                     *uiSessmenu;
  FXMenuPane                     *uiActionmenu;
  FXToolTip                      *m_toolTip;               // controls tooltip behaviour
  FXComboBox                     *m_uiDBSelection;         // Showing current connection
  TableColumnList                *uiTblColList;            // List of tables and columns
  QueryFrame                     *searchFrame;             // Frame to perform queries
  SQLBuffer                      *uiEditor;                // SQL editor
  FXbool                          m_shactive;              // Syntax Highlighting set or not
  LogText                        *uiLog;                   // Result log
  ResultTable                    *uiTable;                 // Result Table
  ConnectionDialog               *m_dgNewCon;              // Connection Dialog
  InsertDialog                   *m_dgInsert;              // Insert Dialog
  
  rdbModel::XercesBuilder        *m_rdbBuilder;            // Builder of the rdb from the xml file
  rdbModel::MysqlConnection      *m_connect;               // Object to connect to a mysql db
  FXString                        m_lastDbSchema;          // last DB schema opened
  rdbModel::Rdb                  *m_rdb;                   // Rdb for last compatible schema
                                                              

  std::string                     m_primaryKey;            // the name of the primary key                       
  
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
    ID_LOWPAN,
    ID_TABLEOUT,
    ID_TOGGLERESULT,
    ID_PASTEFROMTABLE,
    ID_ABOUT,
    ID_INSERT,
    ID_MULTI,
    ID_UPDATELAST,
    ID_UPDATEROW,
    ID_COPYROW,
    ID_INSERTLATEST      
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
  long onUpdResTableCols(FXObject*,FXSelector, void*);
  long onSendQuery(FXObject*,FXSelector, void*);
  long onInsert(FXObject*,FXSelector, void*);
  long onInsertLatest(FXObject*,FXSelector, void*);
  long onMultiInsert(FXObject*,FXSelector, void*);
  long onUpdateLastRow(FXObject*,FXSelector, void*);
  long onUpdateRowByKey(FXObject*,FXSelector, void*);
  long onCopyRowByKey(FXObject*,FXSelector, void* ptr);
  
  // ..eccetera
  
public:

  // RdbGUIWindow constructor
  RdbGUIWindow(FXApp* a);

  // Initialize
  virtual void create();
  
  // RdbGUIWindow destructor
  virtual ~RdbGUIWindow(); 

  // Get the primary key name
  std::string getPrimaryName(){return m_primaryKey;};

  // Set the primary key name
  void setPrimaryName(std::string name){m_primaryKey = name;};
  
private:

  void loadXMLFile(FXString);
  
  void closeConnection();
  };
  
  
#endif

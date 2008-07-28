// $Header$

#include "RdbGuiWindow.h"
#include "InsertDialog.h"

#include "rdbModel/Db/ResultHandle.h"
#include "facilities/commonUtilities.h"

#include <vector>

// Message Map RdbGUIWindow class
FXDEFMAP(RdbGUIWindow) RdbGUIWindowMap[]={

  //__Message_Type____________________ID______________________________Message_Handler_____
  FXMAPFUNC(SEL_UPDATE,   RdbGUIWindow::ID_TITLE,           RdbGUIWindow::onUpdTitle),
  FXMAPFUNC(SEL_UPDATE,   RdbGUIWindow::ID_OVERSTRIKE,      RdbGUIWindow::onUpdOverstrike),
  FXMAPFUNC(SEL_CHANGED,  RdbGUIWindow::ID_SQLEDIT,         RdbGUIWindow::onSQLEdit),
  FXMAPFUNC(SEL_COMMAND,  RdbGUIWindow::ID_OPENXML,         RdbGUIWindow::onOpenXMLFile),
  FXMAPFUNC(SEL_CLOSE,    RdbGUIWindow::ID_TITLE,           RdbGUIWindow::onQuit),
  FXMAPFUNC(SEL_SIGNAL,   RdbGUIWindow::ID_QUIT,            RdbGUIWindow::onQuit),
  FXMAPFUNC(SEL_COMMAND,  RdbGUIWindow::ID_QUIT,            RdbGUIWindow::onQuit),
  FXMAPFUNC(SEL_COMMAND,  RdbGUIWindow::ID_OPENCONNECTION,  RdbGUIWindow::onOpenConnection),
  FXMAPFUNC(SEL_COMMAND,  RdbGUIWindow::ID_CLOSECONNECTION, RdbGUIWindow::onCloseConnection), 
  FXMAPFUNC(SEL_COMMAND,  TableColumnList::ID_TBLLIST,      RdbGUIWindow::onQueryFrameUpdate),
  FXMAPFUNC(SEL_SELECTED, TableColumnList::ID_TBLLIST,      RdbGUIWindow::onUpdResTableCols),
  FXMAPFUNC(SEL_DESELECTED, TableColumnList::ID_TBLLIST,    RdbGUIWindow::onUpdResTableCols),
  FXMAPFUNC(SEL_COMMAND,  QueryFrame::ID_QUERY,             RdbGUIWindow::onSendQuery),
  FXMAPFUNC(SEL_COMMAND,  RdbGUIWindow::ID_INSERT,          RdbGUIWindow::onInsert),
  FXMAPFUNC(SEL_COMMAND,  RdbGUIWindow::ID_MULTI,           RdbGUIWindow::onMultiInsert),  
  FXMAPFUNC(SEL_COMMAND,  RdbGUIWindow::ID_UPDATELAST,      RdbGUIWindow::onUpdateLastRow),
  FXMAPFUNC(SEL_COMMAND,  ResultTable::ID_UPDATEROW,        RdbGUIWindow::onUpdateRowByKey),
  FXMAPFUNC(SEL_COMMAND,  ResultTable::ID_COPYROW,          RdbGUIWindow::onCopyRowByKey),
  FXMAPFUNC(SEL_COMMAND,  RdbGUIWindow::ID_INSERTLATEST,    RdbGUIWindow::onInsertLatest),
  FXMAPFUNC(SEL_COMMAND,  ResultTable::ID_COPYLATEST,    RdbGUIWindow::onCopyLatestByKey),    // added jrb 23 aug 2006

};



// Macro for the GLViewWindow class hierarchy implementation
FXIMPLEMENT(RdbGUIWindow,FXMainWindow,RdbGUIWindowMap,ARRAYNUMBER(RdbGUIWindowMap))


/*******************************************************************************/

// Construct a RdbGUIWindow
  RdbGUIWindow::RdbGUIWindow(FXApp* a):FXMainWindow(a,"rdbGUI",NULL,NULL,DECOR_ALL,0,0,800,600), m_rdb(0)
{
  // Main window set itself as the target
  setTarget(this);
  setSelector(ID_TITLE);

  m_primaryKey = "";

  // Menubar
  uiMenuBar = new FXMenuBar(this, LAYOUT_SIDE_TOP|LAYOUT_FILL_X);

  // Statusbar
  new FXStatusBar(this,
    LAYOUT_SIDE_BOTTOM|LAYOUT_FILL_X|STATUSBAR_WITH_DRAGCORNER);

  // Tooltip
  new FXToolTip(getApp());

  // Content frame
  FXVerticalFrame *uiContent = new FXVerticalFrame(this, LAYOUT_FILL_X|LAYOUT_FILL_Y);

  // Toolbar
  new FXHorizontalSeparator(uiContent, SEPARATOR_RIDGE|LAYOUT_FILL_X);
  FXHorizontalFrame *uiToolbar = new FXHorizontalFrame(uiContent,LAYOUT_FILL_X);

  // DB selection
  new FXLabel(uiToolbar, "Database:");
  m_uiDBSelection = new FXComboBox(uiToolbar, 35, NULL, 0, COMBOBOX_INSERT_LAST|
      FRAME_SUNKEN|FRAME_THICK);

  // Toolbar buttons
  new FXVerticalSeparator(uiToolbar, SEPARATOR_GROOVE|LAYOUT_FILL_Y);
  

  // Horizontal splitter
  FXSplitter *uiHsplitter = new FXSplitter(uiContent, LAYOUT_FILL_X|
                                           LAYOUT_FILL_Y|SPLITTER_TRACKING|SPLITTER_HORIZONTAL);

  // tables and columns frame left
  FXHorizontalFrame *uiTblColframe = new FXHorizontalFrame(uiHsplitter,
                                      LAYOUT_FILL_X|LAYOUT_FILL_Y|FRAME_THICK);
 
  // Treelist
  uiTblColList = new TableColumnList(uiTblColframe, this);

  // Vertical splitter right
  FXSplitter *uiVsplitter = new FXSplitter(uiHsplitter, (LAYOUT_FILL_X|LAYOUT_FILL_Y|
                               SPLITTER_TRACKING|SPLITTER_VERTICAL|PACK_UNIFORM_HEIGHT));


  // Search Frame
  searchFrame = new QueryFrame(uiVsplitter, this);


  // Editor
  //uiEditor = new SQLBuffer(uiEditorframe, this, ID_SQLEDIT, LAYOUT_FILL_X|LAYOUT_FILL_Y);

  FXTabBook *lowerTab = new FXTabBook(uiVsplitter, this, ID_LOWPAN, TABBOOK_BOTTOMTABS|LAYOUT_FILL_X|LAYOUT_FILL_Y);  
  
  new FXTabItem(lowerTab, "Query output", NULL, TAB_BOTTOM_NORMAL);
  
  FXHorizontalFrame *tableFrame = new FXHorizontalFrame(lowerTab,
      LAYOUT_FILL_X|LAYOUT_FILL_Y|FRAME_THICK|FRAME_SUNKEN);
    // Result table
  uiTable = new ResultTable(tableFrame, this, ID_TABLEOUT, LAYOUT_FILL_X|LAYOUT_FILL_Y, 0,0,0,0, 2,2,2,2);
  
  new FXButton(uiToolbar, "Copy", NULL, uiTable, FXTable::ID_COPY_SEL);
      
  new FXTabItem(lowerTab, "Log", NULL, TAB_BOTTOM_NORMAL);         
  // Log window frame
  FXHorizontalFrame *logFrame = new FXHorizontalFrame(lowerTab,
      LAYOUT_FILL_X|LAYOUT_FILL_Y|FRAME_THICK);
  // Result log
  uiLog = new LogText(logFrame, NULL, 0, LAYOUT_FILL_X|LAYOUT_FILL_Y);
  
                   
  // File menu
  uiFilemenu = new FXMenuPane(this);
  new FXMenuTitle(uiMenuBar, "&File", NULL, uiFilemenu);
  new FXMenuCommand(uiFilemenu, "&Open DB Schema ...\tCtl-O\tLoad the xml database description",
    NULL, this, ID_OPENXML);
  new FXMenuSeparator(uiFilemenu);
  new FXMenuCommand(uiFilemenu, "&Quit\tCtl-Q\tQuit DbGui", NULL,
    this, ID_QUIT);
    
  // Session menu
  uiSessmenu = new FXMenuPane(this);
  new FXMenuTitle(uiMenuBar, "&Session", NULL, uiSessmenu);
  m_cmdOpenConn = new FXMenuCommand(uiSessmenu, "&Open connection\tCtl-[\tOpen new database connection...",
      NULL, this, ID_OPENCONNECTION);
  m_cmdCloseConn = new FXMenuCommand(uiSessmenu, "&Close connection\tCtl-]\tClose current database connection",
      NULL, this, ID_CLOSECONNECTION);

  // Action menu
  uiActionmenu = new FXMenuPane(this);
  new FXMenuTitle(uiMenuBar, "&Action", NULL, uiActionmenu);
  new FXMenuCommand(uiActionmenu, "&Insert\tCtl-I\tInsert a new row\tInsert a new row",
    NULL, this, ID_INSERT);
  new FXMenuCommand(uiActionmenu, "&InsertLatest\tCtl-L\tInsert latest row of this type\tInsert latest row of this type",
    NULL, this, ID_INSERTLATEST);
  new FXMenuCommand(uiActionmenu, "&MultiInsert\tCtl-M\tMultiInsert a new row\tMultiInsert a new row",
    NULL, this, ID_MULTI);
  new FXMenuCommand(uiActionmenu, "&Redo Last Insert\tCtl-R\tRedo the last row inserted\tRedo the last row inserted",
    NULL, this, ID_UPDATELAST);

  // Force some reasonable sizes to the layout manager
  uiTblColframe->setWidth(150);
  searchFrame->setHeight(200);
  //uiEditorframe->setHeight(300);
  m_toolTip = new FXToolTip(getApp(), TOOLTIP_VARIABLE);

  // Editor initialization
  //  @filename = 'untitled.sql'
  //  @filenameset = false
  //  @ui_log.editable = false
  //uiEditor->setHiliteMatchTime(1500);
  //uiEditor->setCursorColor(fxcolorfromname("DarkGray"));
  //uiEditor->setFocus();
  //uiEditor->initSyntaxHighlighting();
  //m_shactive = getApp()->reg().readStringEntry("SETTINGS", "HighLightSyntax", "Y") == "Y";
  
  // Result table initialization
  uiTable->setBackColors(FXRGB(255, 255, 255), FXRGB(255, 240, 240),
      FXRGB(240, 255, 240), FXRGB(240, 240, 255));
  uiTable->setTableSize(1, 1);
  uiTable->setRowHeaderWidth(20);
  uiTable->setItemText(0, 0, "No data");
  
  // Disable some menù items
  m_cmdOpenConn->disable();
  m_cmdCloseConn->disable();
  
  // Initialize connection dialog and mysql connection
  m_dgNewCon = new ConnectionDialog(this);
  m_connect = new rdbModel::MysqlConnection(uiLog->getOutStream(), uiLog->getErrStream());
  
  // Initialize insert dialog 
  m_dgInsert = new InsertDialog(getApp());

  
  // Initialize the rdb builder
  m_rdbBuilder = new rdbModel::XercesBuilder();
  
  m_lastDbSchema = "";
}


// Destructor
RdbGUIWindow::~RdbGUIWindow()
{
  delete uiMenuBar;
  delete uiFilemenu;
  delete uiSessmenu;
  delete uiActionmenu;
  delete m_uiDBSelection;
}



// Create and initialize
void RdbGUIWindow::create()
{
  FXMainWindow::create();
  show(PLACEMENT_SCREEN);
}


// Update title
long  RdbGUIWindow::onUpdTitle(FXObject*, FXSelector , void*)
{
//   FXString title = "RdbGUI - [" + filename + "]";
//   if (uiEditor.isModified())
//     title << "*";
//   sender.handle(this, MKUINT(FXWindow::ID_SETSTRINGVALUE, SEL_COMMAND), title);
  return 1;
}

// Update box for overstrike mode display
long  RdbGUIWindow::onUpdOverstrike(FXObject*, FXSelector, void*)
{
//   mode = ((@ui_editor.textStyle & TEXT_OVERSTRIKE) != 0) ? "OVR" : "INS";
//   sender.handle(self, MKUINT(ID_SETSTRINGVALUE, SEL_COMMAND), mode);
  return 1;
}

long RdbGUIWindow::onSQLEdit(FXObject*, FXSelector, void*)
{
  if (m_shactive)
    uiEditor->highlightSyntax();
  return 1;
}

long RdbGUIWindow::onOpenXMLFile(FXObject*, FXSelector, void*)
{
  // Open file
  FXFileDialog *opendialog = new FXFileDialog(this, "Open File");
  opendialog->setSelectMode(SELECTFILE_EXISTING);
  opendialog->setPatternList("XML files (*.xml)\nAll files (*)"); 

  // If the env variable RDBMODELROOT is set, the directory $RDBMODELROOT/xml is
  // used as the default opening directory for the file browser
  if (facilities::commonUtilities::getXmlPath("rdbModel") != "") {
    std::string rdbModelXml = 
      facilities::commonUtilities::getXmlPath("rdbModel");
    opendialog->setDirectory(FXString(rdbModelXml.c_str()));
  }
  
  if (opendialog->execute() != 0)
    {
      if (onCloseConnection(NULL,0,NULL))
        {
          if (m_rdb) {
            delete m_rdb;
          }
          m_rdb = new rdbModel::Rdb;

          loadXMLFile(opendialog->getFilename());        
        }
    }
    
  return 1;
}


void RdbGUIWindow::loadXMLFile(FXString fileName)
{
  if (FXFile::exists(fileName))
    {
      if (m_rdb->build(fileName.text(), m_rdbBuilder)) 
        {
          FXMessageBox::error(this, MBOX_OK, "Error: rdbModel building", 
              "an error as occurred during the construction of the rdb model");
          m_cmdOpenConn->disable();
        }
      else
        {
          m_rdb->accept(uiTblColList);
          m_cmdOpenConn->enable();
          m_lastDbSchema = fileName;
        }
        
    } 
  searchFrame->setEnabled(false); 
}


long RdbGUIWindow::onQuit(FXObject*, FXSelector, void*)
{
  getApp()->exit(0);
  return 1;
}

long RdbGUIWindow::onOpenConnection(FXObject*,FXSelector, void*)
{

  if (m_dgNewCon->execute(PLACEMENT_OWNER) != 0)
    {
      getApp()->beginWaitCursor();
      std::vector<FXString> data = m_dgNewCon->getConnectionData();      
      if (!(m_connect->open(data[1].text(), data[2].text(), data[3].text(), data[0].text()))) 
        {
          uiLog->update();
        }
      else
        {
          uiLog->update();
          m_uiDBSelection->appendItem(data[0]+" ("+data[2]+"@"+data[1]+")");
          m_uiDBSelection->setCurrentItem(0);   // this line will have to be changed
          
          rdbModel::MATCH match = m_connect->matchSchema(m_rdb, false);
        
          switch (match) {
          case rdbModel::MATCHequivalent:
            uiLog->update();   
            uiLog->logText("XML schema and MySQL database are equivalent!\n");
            searchFrame->setConnection(m_connect);
            m_cmdOpenConn->disable();
            m_cmdCloseConn->enable();
            break;
          case rdbModel::MATCHcompatible:
            uiLog->update();   
            uiLog->logText("XML schema and MySQL database are compatible\n");
            searchFrame->setConnection(m_connect);
            m_cmdOpenConn->disable();
            m_cmdCloseConn->enable();
            break;
          case rdbModel::MATCHfail:
            uiLog->update();   
            uiLog->logText("XML schema and MySQL database are NOT compatible\n");

            if (m_rdb) {
              delete m_rdb;
              m_rdb = 0;
            }
            m_connect->close();
            break;
          case rdbModel::MATCHnoConnection:
            uiLog->update();   
            uiLog->logText("Connection failed while attempting match\n");
            if (m_rdb) {
              delete m_rdb;
              m_rdb = 0;
            }
            m_connect->close();
            break;
          }
          
        }
      getApp()->endWaitCursor();
    }
  return 1;
}


long RdbGUIWindow::onCloseConnection(FXObject*,FXSelector, void*)
{
  if (m_uiDBSelection->getNumItems() == 0)
    return 1;

  FXString message = "Close connection to " + m_uiDBSelection->getText()+"?\n";
  if (MBOX_CLICKED_OK == FXMessageBox::question(this, MBOX_OK_CANCEL, "Close Connection",
      message.text()))
    {
      closeConnection();
      m_cmdCloseConn->disable();
      m_cmdOpenConn->enable();
      return 1;
    }
  return 0;  
}

long RdbGUIWindow::onUpdResTableCols(FXObject*,FXSelector, void* ptr)
{
  const FXCheckList* chList = uiTblColList->getColList();
  FXint index = (FXint) ptr;
  
  if (index >= 0 && index < chList->getNumItems());
    {
 
      
      if (!chList->isItemChecked(index))
        uiTable->hideColumn(index);
      else
        uiTable->showColumn(index);
    }
  return 1;
}

long RdbGUIWindow::onSendQuery(FXObject*,FXSelector, void*)
{
  FXCheckList *columns = (FXCheckList *) uiTblColList->getColList();
  rdbModel::ResultHandle *queryResult = searchFrame->getQueryResult();
  uiLog->update();

  int index = uiTblColList->getTableList()->getCurrentItem(); 
  uiTable->setTableName((uiTblColList->getTableList()->getItemText(index)).text());

  if (!queryResult) 

    {
      uiTable->setTableSize(1, 1);
      uiTable->setRowHeaderWidth(20);
      uiTable->setItemText(0, 0, "Null query result");
      return 1;
    }

  if (queryResult->getNRows() < 1)
    {
      uiTable->setTableSize(1, 1);
      uiTable->setRowHeaderWidth(20);
      uiTable->setItemText(0, 0, "No data");
      return 1;
    }
    
  uiTable->setTableSize(queryResult->getNRows(), columns->getNumItems());
  

  int i,j;
  for (j = 0; j < uiTable->getNumColumns(); j++)
    {
      uiTable->setColumnText(j, columns->getItemText(j));
      uiTable->getColumnHeader()->getItem(j)->setIcon(columns->getItemIcon(j));
    }
    
  for (i = 0; i < uiTable->getNumRows(); i++)
    uiTable->setRowText(i, FXStringVal(i+1));


  std::vector<std::string> rowValues;
  for (i=0; i < uiTable->getNumRows(); i++)
    {
      queryResult->getRow(rowValues, i, 1);
      for (j = 0; j < uiTable->getNumColumns(); j++)
        {
          uiTable->setItemText(i, j, rowValues[j].c_str());
        }
    }
  uiTable->format();
  
  for (i = 0; i < columns->getNumItems(); i++)
    onUpdResTableCols(NULL, 0, (void *) i);
  return 1;
}

long RdbGUIWindow::onMultiInsert(FXObject*,FXSelector, void*)
{
  m_dgInsert->setConnection(m_connect);
  m_dgInsert->setRdb(m_rdb);
  if ((uiTblColList->getTableList()->getNumItems() == 0 )
      || (m_connect == 0) || !(m_connect->isConnected()) ||
      (m_dgInsert->getLastRow() < 0))
    return 1;

  int index = uiTblColList->getTableList()->getCurrentItem(); 
  m_dgInsert->setTableName((uiTblColList->getTableList()->getItemText(index)).text());
  // Fill the form
  m_rdb->accept(m_dgInsert);
  // Fill the form with the sticky values of the last row inserted in the DB
  m_dgInsert->fillStickyWithLastRow();
  // Set the mode of the dialog to Insert
  m_dgInsert->setInsertMode(1);
  
  m_dgInsert->show(PLACEMENT_OWNER);
  m_dgInsert->setUiLog(uiLog);

  if (m_dgInsert->getDefaultWidth() < 250)
    m_dgInsert->resize(250,m_dgInsert->getDefaultHeight());
  else
    m_dgInsert->resize(m_dgInsert->getDefaultWidth(),m_dgInsert->getDefaultHeight());
  m_dgInsert->recalc();  
  
  return 1;
}


long RdbGUIWindow::onInsert(FXObject*,FXSelector, void*)
{
  m_dgInsert->setConnection(m_connect);
  m_dgInsert->setRdb(m_rdb);
  if ((uiTblColList->getTableList()->getNumItems() == 0 )
      || (m_connect == 0) || !(m_connect->isConnected()))
    return 1;

  int index = uiTblColList->getTableList()->getCurrentItem(); 
  m_dgInsert->setTableName((uiTblColList->getTableList()->getItemText(index)).text());
  // Fill the form
  m_rdb->accept(m_dgInsert);

  // Set the mode of the dialog to Insert
  m_dgInsert->setInsertMode(1);
  
  m_dgInsert->show(PLACEMENT_OWNER);
  m_dgInsert->setUiLog(uiLog);

  if (m_dgInsert->getDefaultWidth() < 250)
    m_dgInsert->resize(250,m_dgInsert->getDefaultHeight());
  else
    m_dgInsert->resize(m_dgInsert->getDefaultWidth(),m_dgInsert->getDefaultHeight());
  m_dgInsert->recalc();  
  
  return 1;
}

long RdbGUIWindow::onInsertLatest(FXObject*,FXSelector, void*)
{
  m_dgInsert->setConnection(m_connect);
  m_dgInsert->setRdb(m_rdb);
  if ((uiTblColList->getTableList()->getNumItems() == 0 )
      || (m_connect == 0) || !(m_connect->isConnected()))
    return 1;

  int index = uiTblColList->getTableList()->getCurrentItem(); 
  m_dgInsert->setTableName((uiTblColList->getTableList()->getItemText(index)).text());
  // Fill the form
  m_rdb->accept(m_dgInsert);
  // Set the mode of the dialog to InsertLatest
  m_dgInsert->setInsertMode(2);
  
  m_dgInsert->show(PLACEMENT_OWNER);
  m_dgInsert->setUiLog(uiLog);

  if (m_dgInsert->getDefaultWidth() < 250)
    m_dgInsert->resize(250,m_dgInsert->getDefaultHeight());
  else
    m_dgInsert->resize(m_dgInsert->getDefaultWidth(),m_dgInsert->getDefaultHeight());
  m_dgInsert->recalc();  
  
  return 1;
}


long RdbGUIWindow::onUpdateLastRow(FXObject*,FXSelector, void*)
{
  m_dgInsert->setConnection(m_connect);
  m_dgInsert->setRdb(m_rdb);
  if ((uiTblColList->getTableList()->getNumItems() == 0 )
      || (m_connect == 0) || !(m_connect->isConnected()) 
      || (m_dgInsert->getLastRow() < 0))
    return 1;

  m_dgInsert->setTableName(m_dgInsert->getLastTblName());
  // Build the form
  m_rdb->accept(m_dgInsert);
  // Fill the form with the last row inserted in the DB
  m_dgInsert->fillWithLastRow();
  // Set the mode of the dialog to Update Last Row
  m_dgInsert->setInsertMode(0);
  
  m_dgInsert->show(PLACEMENT_OWNER);
  m_dgInsert->setUiLog(uiLog);

  if (m_dgInsert->getDefaultWidth() < 250)
    m_dgInsert->resize(250,m_dgInsert->getDefaultHeight());
  else
    m_dgInsert->resize(m_dgInsert->getDefaultWidth(),m_dgInsert->getDefaultHeight());
  m_dgInsert->recalc();  
  
  return 1;
}

long RdbGUIWindow::onUpdateRowByKey(FXObject*,FXSelector, void* ptr)
{
  FXint row = (FXint) ptr;
  
  m_dgInsert->setConnection(m_connect);
  m_dgInsert->setRdb(m_rdb);
  if ((uiTblColList->getTableList()->getNumItems() == 0 )
      || (m_connect == 0) || !(m_connect->isConnected()))
    return 1;

  m_dgInsert->setTableName(uiTable->getTableName());
  // Build the form
  m_rdb->accept(m_dgInsert);
  // Fill the form with the last row inserted in the DB
  m_dgInsert->fillWithRowByKey(uiTable->getItemText(row, 0).text());
  // Set the mode of the dialog to Update Last Row
  m_dgInsert->setInsertMode(0);
  
  m_dgInsert->show(PLACEMENT_OWNER);
  m_dgInsert->setUiLog(uiLog);

  if (m_dgInsert->getDefaultWidth() < 250)
    m_dgInsert->resize(250,m_dgInsert->getDefaultHeight());
  else
    m_dgInsert->resize(m_dgInsert->getDefaultWidth(),m_dgInsert->getDefaultHeight());
  m_dgInsert->recalc();  
  
  return 1;
}

long RdbGUIWindow::onCopyRowByKey(FXObject*,FXSelector, void* ptr)
{
  FXint row = (FXint) ptr;
  
  m_dgInsert->setConnection(m_connect);
  m_dgInsert->setRdb(m_rdb);
  if ((uiTblColList->getTableList()->getNumItems() == 0 )
      || (m_connect == 0) || !(m_connect->isConnected()))
    return 1;

  m_dgInsert->setTableName(uiTable->getTableName());
  // Build the form
  m_rdb->accept(m_dgInsert);
  // Fill the form with the last row inserted in the DB
  m_dgInsert->fillWithRowByKey(uiTable->getItemText(row, 0).text());
  // Set the mode of the dialog to Insert
  m_dgInsert->setInsertMode(1);
  
  m_dgInsert->show(PLACEMENT_OWNER);
  m_dgInsert->setUiLog(uiLog);

  if (m_dgInsert->getDefaultWidth() < 250)
    m_dgInsert->resize(250,m_dgInsert->getDefaultHeight());
  else
    m_dgInsert->resize(m_dgInsert->getDefaultWidth(),m_dgInsert->getDefaultHeight());
  m_dgInsert->recalc();  
  
  return 1;
}

 // Insert a row similar to selected one, but invoke rdbModel::insertLatest rather than just insert
long RdbGUIWindow::onCopyLatestByKey(FXObject*,FXSelector, void* ptr)
{
  FXint row = (FXint) ptr;
  
  m_dgInsert->setConnection(m_connect);
  m_dgInsert->setRdb(m_rdb);
  if ((uiTblColList->getTableList()->getNumItems() == 0 )
      || (m_connect == 0) || !(m_connect->isConnected()))
    return 1;

  m_dgInsert->setTableName(uiTable->getTableName());
  // Build the form
  m_rdb->accept(m_dgInsert);
  // Fill the form with the last row inserted in the DB
  m_dgInsert->fillWithRowByKey(uiTable->getItemText(row, 0).text());
  // Set the mode of the dialog to InsertLatest 
  m_dgInsert->setInsertMode(2);
  
  m_dgInsert->show(PLACEMENT_OWNER);
  m_dgInsert->setUiLog(uiLog);

  if (m_dgInsert->getDefaultWidth() < 250)
    m_dgInsert->resize(250,m_dgInsert->getDefaultHeight());
  else
    m_dgInsert->resize(m_dgInsert->getDefaultWidth(),m_dgInsert->getDefaultHeight());
  m_dgInsert->recalc();  
  
  return 1;
}

void RdbGUIWindow::closeConnection()
{
  m_connect->close();

  m_uiDBSelection->removeItem(m_uiDBSelection->getCurrentItem());
  searchFrame->reset();
  uiTable->clearItems();
  uiTable->setTableSize(1, 1);
  uiTable->setItemText(0, 0, "No data");
  uiTable->setTableName("");
}

long RdbGUIWindow::onQueryFrameUpdate(FXObject *, FXSelector, void*)
{
  searchFrame->setEnabled(true);
  searchFrame->updateColumnSelection(uiTblColList->getTableList(), 
      uiTblColList->getColList());
  return 1; 
}







#include "RdbGuiWindow.h"
#include "InsertDialog.h"

#include "rdbModel/Db/ResultHandle.h"

#include <vector>

// Message Map RdbGUIWindow class
FXDEFMAP(RdbGUIWindow) RdbGUIWindowMap[]={

  //__Message_Type_____________ID________________________Message_Handler_____
  FXMAPFUNC(SEL_UPDATE,   RdbGUIWindow::ID_TITLE,              RdbGUIWindow::onUpdTitle),
  FXMAPFUNC(SEL_UPDATE,   RdbGUIWindow::ID_OVERSTRIKE,         RdbGUIWindow::onUpdOverstrike),
  FXMAPFUNC(SEL_CHANGED,  RdbGUIWindow::ID_SQLEDIT,            RdbGUIWindow::onSQLEdit),
  FXMAPFUNC(SEL_COMMAND,  RdbGUIWindow::ID_OPENXML,            RdbGUIWindow::onOpenXMLFile),
  FXMAPFUNC(SEL_CLOSE,    RdbGUIWindow::ID_TITLE,              RdbGUIWindow::onQuit),
  FXMAPFUNC(SEL_SIGNAL,   RdbGUIWindow::ID_QUIT,               RdbGUIWindow::onQuit),
  FXMAPFUNC(SEL_COMMAND,  RdbGUIWindow::ID_QUIT,               RdbGUIWindow::onQuit),
  FXMAPFUNC(SEL_COMMAND,  RdbGUIWindow::ID_OPENCONNECTION,     RdbGUIWindow::onOpenConnection),
  FXMAPFUNC(SEL_COMMAND,  RdbGUIWindow::ID_CLOSECONNECTION,    RdbGUIWindow::onCloseConnection), 
  FXMAPFUNC(SEL_COMMAND,  TableColumnList::ID_TBLLIST,         RdbGUIWindow::onQueryFrameUpdate), 
  FXMAPFUNC(SEL_COMMAND,  QueryFrame::ID_QUERY,                RdbGUIWindow::onSendQuery),
  FXMAPFUNC(SEL_COMMAND,  RdbGUIWindow::ID_INSERT,             RdbGUIWindow::onInsert)
};



// Macro for the GLViewWindow class hierarchy implementation
FXIMPLEMENT(RdbGUIWindow,FXMainWindow,RdbGUIWindowMap,ARRAYNUMBER(RdbGUIWindowMap))


/*******************************************************************************/

// Construct a RdbGUIWindow
RdbGUIWindow::RdbGUIWindow(FXApp* a):FXMainWindow(a,"rdbGUI",NULL,NULL,DECOR_ALL,0,0,800,600)
{

  // Main window set itself as the target
  setTarget(this);
  setSelector(ID_TITLE);

  // Menubar
  uiMenuBar = new FXMenuBar(this, LAYOUT_SIDE_TOP|LAYOUT_FILL_X);

  // Statusbar
  FXStatusBar *uiStatusBar = new FXStatusBar(this,
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
  new FXToggleButton(uiToolbar, "&Tab", "&Txt\tToggle output",
                     NULL, NULL, this, ID_TOGGLERESULT, FRAME_RAISED|FRAME_THICK);
  new FXButton(uiToolbar, "&Paste result\tPaste from result table",
               NULL, this, ID_PASTEFROMTABLE, FRAME_RAISED|FRAME_THICK);
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
      
  new FXTabItem(lowerTab, "Log", NULL, TAB_BOTTOM_NORMAL);         
  // Log window frame
  FXHorizontalFrame *logFrame = new FXHorizontalFrame(lowerTab,
      LAYOUT_FILL_X|LAYOUT_FILL_Y|FRAME_THICK);
  // Result log
  uiLog = new LogText(logFrame, NULL, 0, LAYOUT_FILL_X|LAYOUT_FILL_Y);
  
                   
  // File menu
  uiFilemenu = new FXMenuPane(this);
  new FXMenuTitle(uiMenuBar, "&File", NULL, uiFilemenu);
  new FXMenuCommand(uiFilemenu, "&Open XML ...\tCtl-O\tLoad the xml database description",
    NULL, this, ID_OPENXML);
  new FXMenuSeparator(uiFilemenu);
  new FXMenuCommand(uiFilemenu, "&Quit\tCtl-Q\tQuit DbGui", NULL,
    this, ID_QUIT);
    
  // Session menu
  uiSessmenu = new FXMenuPane(this);
  new FXMenuTitle(uiMenuBar, "&Session", NULL, uiSessmenu);
  new FXMenuCommand(uiSessmenu, "&Open connection\tCtl-[\tOpen new database connection...",
      NULL, this, ID_OPENCONNECTION);
  new FXMenuCommand(uiSessmenu, "&Close connection\tCtl-]\tClose current database connection",
      NULL, this, ID_CLOSECONNECTION);

  // Action menu
  uiActionmenu = new FXMenuPane(this);
  new FXMenuTitle(uiMenuBar, "&Action", NULL, uiActionmenu);
  new FXMenuCommand(uiActionmenu, "&Insert\tCtl-I\tInsert a new row\tInsert a new row",
    NULL, this, ID_INSERT);

  // Force some reasonable sizes to the layout manager
  uiTblColframe->setWidth(150);
  searchFrame->setHeight(200);
  //uiEditorframe->setHeight(300);

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
  uiTable->setTableSize(1, 1);
  uiTable->setItemText(0, 0, "No data");
  
  
  // Initialize connection dialog and mysql connection
  m_dgNewCon = new ConnectionDialog(this);
  m_connect = new rdbModel::MysqlConnection(uiLog->getOutStream(), uiLog->getErrStream());

  // Initialize insert dialog 
  m_dgInsert = new InsertDialog(this);

  
  // Initialize the rdb manager and it's builder
  m_rdbBuilder = new rdbModel::XercesBuilder();
  m_rdbManager = rdbModel::Manager::getManager();
  m_rdbManager->setBuilder(m_rdbBuilder);

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
  onCloseConnection(NULL,0,NULL);
  if (m_rdbManager->getRdb())
    m_rdbManager->cleanRdb();
  FXFileDialog *opendialog = new FXFileDialog(this, "Open File");
  opendialog->setSelectMode(SELECTFILE_EXISTING);
  opendialog->setPatternList("XML files (*.xml)\nAll files (*)");
  if (opendialog->execute() != 0)
    loadXMLFile(opendialog->getFilename());
  return 1;
}


void RdbGUIWindow::loadXMLFile(FXString fileName)
{
  if (FXFile::exists(fileName))
    {
      m_rdbManager->setInputSource(fileName.text());
      if (m_rdbManager->build())
        FXMessageBox::error(this, MBOX_OK, "Error: rdbModel building", 
                            "an error as occurred during the construction of the rdb model");
      m_rdbManager->startVisitor(uiTblColList);
    }  
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
          
          rdbModel::MATCH match = m_connect->matchSchema(m_rdbManager->getRdb());
        
          switch (match) {
          case rdbModel::MATCHequivalent:
            uiLog->logText("XML schema and MySQL database are equivalent!\n");
            searchFrame->setConnection(m_connect);
            break;
          case rdbModel::MATCHcompatible:
            uiLog->logText("XML schema and MySQL database are compatible\n");
            searchFrame->setConnection(m_connect);
            break;
          case rdbModel::MATCHfail:
            uiLog->logText("XML schema and MySQL database are NOT compatible\n");
            break;
          case rdbModel::MATCHnoConnection:
            uiLog->logText("Connection failed while attempting match\n");
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
      m_connect->close();
      m_uiDBSelection->removeItem(m_uiDBSelection->getCurrentItem());
    }
  return 1;  
}

long RdbGUIWindow::onSendQuery(FXObject*,FXSelector, void*)
{
  FXCheckList *columns = (FXCheckList *) uiTblColList->getColList();
  rdbModel::ResultHandle *queryResult = searchFrame->getQueryResult();
  uiTable->setTableSize(queryResult->getNRows(), columns->getNumItems());
  
  int i,j;
  for (j = 0; j < uiTable->getNumColumns(); j++)
    uiTable->setColumnText(j, columns->getItemText(j));
    
  for (i = 0; i < uiTable->getNumRows(); i++)
    uiTable->setRowText(i, FXStringVal(i+1));
    
  for (i=0; i < uiTable->getNumRows(); i++)
    {
      std::vector<std::string> rowValues;
      queryResult->getRow(rowValues, i);
      for (j = 0; j < uiTable->getNumColumns(); j++)
        {
          uiTable->setItemText(i, j, rowValues[j].c_str());
        }
    }
  return 1;
}

long RdbGUIWindow::onInsert(FXObject*,FXSelector, void*)
{
  m_dgInsert->setConnection(m_connect);
  m_dgInsert->setTableName("metadata_v2r1");
  
  m_rdbManager->startVisitor(m_dgInsert);
  m_dgInsert->show();
  m_dgInsert->update();

  if (m_dgInsert->execute(PLACEMENT_OWNER) != 0)
  {
    return 1; 
  }

  uiLog->update();
  
  return 0;
}

void RdbGUIWindow::closeConnection()
{
  m_connect->close();
  m_uiDBSelection->removeItem(m_uiDBSelection->getCurrentItem());
  
}

long RdbGUIWindow::onQueryFrameUpdate(FXObject *, FXSelector, void*)
{
  searchFrame->updateColumnSelection(uiTblColList->getTableList(), 
      uiTblColList->getColList());
  return 1; 
}


// Here we begin
int main(int argc,char *argv[])
{

  // Make application
  FXApp application("RdbGUI","Calibration");

  // Open the display
  application.init(argc,argv);

  // Make window
  new RdbGUIWindow(&application);

  // Create the application's windows
  application.create();

  // Run the application
  return application.run();
}





#include "TableColumnList.h"
#include "Icons.h"
#include "RdbGuiWindow.h"

// Message Map TableColumnList class
FXDEFMAP(TableColumnList) TableColumnListMap[]={

  //__Message_Type_____________________________ID________________________Message_Handler_____
  FXMAPFUNC(SEL_COMMAND,           TableColumnList::ID_TBLLIST,          TableColumnList::onSelectTable),
  FXMAPFUNC(SEL_SELECTED,          TableColumnList::ID_COLLIST,          TableColumnList::onCheckColumn),
  FXMAPFUNC(SEL_DESELECTED,        TableColumnList::ID_COLLIST,          TableColumnList::onUncheckColumn),
  FXMAPFUNC(SEL_RIGHTBUTTONRELEASE,TableColumnList::ID_COLLIST,          TableColumnList::onCmdMenuPane)
  };

// Object implementation
FXIMPLEMENT(TableColumnList,FXVerticalFrame,TableColumnListMap,ARRAYNUMBER(TableColumnListMap))


TableColumnList::TableColumnList(FXComposite *owner, RdbGUIWindow *target, FXSelector sel):
  FXVerticalFrame(owner, LAYOUT_FILL_X|LAYOUT_FILL_Y|FRAME_SUNKEN, 0, 0, 0, 0, 0, 0, 0, 0), 
  m_target(target), 
  m_selector(sel)
{

 // Vertical splitter right
  FXSplitter *uiVsplitter = new FXSplitter(this, (LAYOUT_FILL_X|LAYOUT_FILL_Y|
      SPLITTER_TRACKING|SPLITTER_VERTICAL));
           
  FXVerticalFrame *tblFrame = new FXVerticalFrame(uiVsplitter,
    LAYOUT_FILL_X|LAYOUT_FILL_Y, 0, 0, 0, 0, 0, 0, 0, 0);
                                           
  new FXLabel(tblFrame, "Tables", new FXICOIcon(getApp(), tablesImg), LABEL_NORMAL,0,0,0,0,4,0,0,0);
  
  FXVerticalFrame *tblListFrame = new FXVerticalFrame(tblFrame, LAYOUT_FILL_X|LAYOUT_FILL_Y|FRAME_SUNKEN, 0, 0, 0, 0, 0, 0, 0, 0);  
  
  m_tblList = new FXList(tblListFrame, this, ID_TBLLIST, LIST_SINGLESELECT|LAYOUT_FILL_X|LAYOUT_FILL_Y);
  
  FXVerticalFrame *colFrame = new FXVerticalFrame(uiVsplitter,
    LAYOUT_FILL_X|LAYOUT_FILL_Y, 0, 0, 0, 0, 0, 0, 0, 0);
                                           
  new FXLabel(colFrame, "Columns", new FXICOIcon(getApp(), columnsImg), LABEL_NORMAL,0,0,0,0,4,0,0,0);
  
  FXVerticalFrame *colListFrame = new FXVerticalFrame(colFrame, LAYOUT_FILL_X|LAYOUT_FILL_Y|FRAME_SUNKEN, 0, 0, 0, 0, 0, 0, 0, 0); 
        
  m_colList = new FXCheckList(colListFrame, 0,this, ID_COLLIST, LIST_NORMAL|LAYOUT_FILL_X|LAYOUT_FILL_Y);
  
  m_colListPop = new FXMenuPane(owner);
  new FXMenuCommand(m_colListPop,"check selected",NULL,m_colList,FXCheckList::ID_CHECKSEL);
  new FXMenuCommand(m_colListPop,"uncheck selected",NULL,m_colList,FXCheckList::ID_UNCHECKSEL); 
  new FXMenuCommand(m_colListPop,"check all",NULL,m_colList,FXCheckList::ID_CHECKALL); 
  new FXMenuCommand(m_colListPop,"uncheck all",NULL,m_colList,FXCheckList::ID_UNCHECKALL); 
      
  m_tableSelected = false;
  
  m_primKeyIcon = new FXGIFIcon(getApp(),primkey);
  m_primKeyIcon->create();
  
//   m_checkIcon = new FXICOIcon(getApp(),checkedImg);
//   m_checkIcon->create();
//   m_uncheckIcon = new FXICOIcon(getApp(),uncheckedImg);
//   m_uncheckIcon->create();
}

long TableColumnList::onSelectTable(FXObject*,FXSelector,void*)
{
  m_tableSelected = true;
  m_colList->clearItems();
  m_target->setPrimaryName("");
  rdbModel::Table *table = (rdbModel::Table*)(m_tblList->getItemData(m_tblList->getCurrentItem()));
  table->accept(this);
  m_tableSelected = false;
  return m_target && m_target->handle(NULL,FXSEL(SEL_COMMAND,ID_TBLLIST),NULL);
}


long TableColumnList::onCheckColumn(FXObject*,FXSelector,void* ptr)
{
  return m_target && m_target->handle(NULL,FXSEL(SEL_SELECTED,ID_TBLLIST),ptr);;
}

long TableColumnList::onUncheckColumn(FXObject*,FXSelector,void* ptr)
{
  return m_target && m_target->handle(NULL,FXSEL(SEL_DESELECTED,ID_TBLLIST),ptr);;
}

long TableColumnList::onCmdMenuPane(FXObject*, FXSelector, void* ptr)
{
  FXEvent *event = (FXEvent*) ptr;
  m_colListPop->popup(NULL,event->root_x,event->root_y);
  //getApp()->runModalWhileShown(getOwner());
  return 1;
}

void TableColumnList::reset()
{
  m_tblList->clearItems();
  m_colList->clearItems(); 
}

rdbModel::Visitor::VisitorState TableColumnList::visitRdb(rdbModel::Rdb *)
{
  m_tblList->clearItems();
  m_colList->clearItems();
  return rdbModel::Visitor::VCONTINUE;
}


rdbModel::Visitor::VisitorState TableColumnList::visitTable(rdbModel::Table *table)
{
  if (m_tableSelected)
    {
      return rdbModel::Visitor::VCONTINUE;
    }
  FXListItem *item = new FXListItem(table->getName().c_str(), NULL, table);
  m_tblList->appendItem(item);
  return rdbModel::Visitor::VBRANCHDONE;
}


rdbModel::Visitor::VisitorState TableColumnList::visitColumn(rdbModel::Column *column)
{
  FXCheckListItem *item = new FXCheckListItem(column->getName().c_str(), NULL, column);
  if (column->isPrimaryKey())
  {
    item->setIcon(m_primKeyIcon);    
    m_target->setPrimaryName(column->getName());
  }
  item->setChecked(true);
  item->setTipText(column->getComment().c_str());
  m_colList->appendItem(item);
  return rdbModel::Visitor::VBRANCHDONE;
}


rdbModel::Visitor::VisitorState TableColumnList::visitIndex(rdbModel::Index *)
{
  return rdbModel::Visitor::VCONTINUE;
}

rdbModel::Visitor::VisitorState TableColumnList::visitAssertion(rdbModel::Assertion *)
{
  return rdbModel::Visitor::VCONTINUE;
}


#include "TableColumnList.h"
#include "Icons.h"

// Message Map TableColumnList class
FXDEFMAP(TableColumnList) TableColumnListMap[]={

  //__Message_Type_____________ID________________________Message_Handler_____
  FXMAPFUNC(SEL_COMMAND,   TableColumnList::ID_TBLLIST,          TableColumnList::onSelectTable),
  FXMAPFUNC(SEL_CLICKED,   TableColumnList::ID_COLLIST,          TableColumnList::onSelectColumn)
  };

// Object implementation
FXIMPLEMENT(TableColumnList,FXVerticalFrame,TableColumnListMap,ARRAYNUMBER(TableColumnListMap))


TableColumnList::TableColumnList(FXComposite *owner, FXObject *target, FXSelector sel):
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
      
  m_tableSelected = false;
  
//   m_checkIcon = new FXICOIcon(getApp(),checkedImg);
//   m_checkIcon->create();
//   m_uncheckIcon = new FXICOIcon(getApp(),uncheckedImg);
//   m_uncheckIcon->create();
}

long TableColumnList::onSelectTable(FXObject*,FXSelector sel,void*)
{
  m_tableSelected = true;
  m_colList->clearItems();
  rdbModel::Table *table = (rdbModel::Table*)(m_tblList->getItemData(m_tblList->getCurrentItem()));
  table->accept(this);
  m_tableSelected = false;
  return m_target && m_target->handle(NULL,FXSEL(SEL_COMMAND,ID_TBLLIST),NULL);
}


long TableColumnList::onSelectColumn(FXObject*,FXSelector,void*)
{
  return 1;
}

rdbModel::Visitor::VisitorState TableColumnList::visitRdb(rdbModel::Rdb *rdb)
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
  item->setChecked(true);
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

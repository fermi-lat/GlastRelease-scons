
#include "QueryFrame.h"
#include "FXCheckList.h"
#include "ColWidgetFactory.h"

// Message Map TableColumnList class
FXDEFMAP(QueryFrame) QueryFrameMap[]={

  //__Message_Type_____________ID________________________Message_Handler_____
  FXMAPFUNC(SEL_COMMAND,   QueryFrame::ID_MORE,          QueryFrame::onCmdMore),
  FXMAPFUNC(SEL_COMMAND,   QueryFrame::ID_FEWER,         QueryFrame::onCmdFewer),
  FXMAPFUNC(SEL_COMMAND,   QueryFrame::ID_COLSELECT,     QueryFrame::onSelectCol)
  
  };

// Object implementation
FXIMPLEMENT(QueryFrame,FXVerticalFrame,QueryFrameMap,ARRAYNUMBER(QueryFrameMap))



QueryFrame::QueryFrame(FXComposite *owner):
  FXVerticalFrame(owner, LAYOUT_FILL_X|LAYOUT_FILL_Y, 0, 0, 0, 0, 0, 0, 0, 0)
{

  // Search frame scroll window
  FXScrollWindow *srcScrollWindow = new FXScrollWindow(this, LAYOUT_FILL_X|LAYOUT_FILL_Y);
  
  // Matrix of FXComboBox 
  m_searchFrame = new FXMatrix(srcScrollWindow, 4, MATRIX_BY_COLUMNS|LAYOUT_FILL_X|
      LAYOUT_FILL_Y|FRAME_SUNKEN|FRAME_THICK|PACK_UNIFORM_WIDTH);
  m_searchFrame->setBackColor(getApp()->getBackColor());
  
  FXComboBox *temp;
  new FXComboBox(m_searchFrame, 0, this, ID_COLSELECT, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_STATIC|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
  new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_STATIC|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
  FXHorizontalFrame *tempFrame = new FXHorizontalFrame(m_searchFrame,LAYOUT_FILL_X|LAYOUT_FILL_Y|LAYOUT_FILL_COLUMN,  
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  tempFrame->setBackColor(getApp()->getBackColor());
  new FXComboBox(tempFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_NORMAL|LAYOUT_FILL_X);
  
  FXHorizontalFrame *addRemovFrame = new FXHorizontalFrame(this, LAYOUT_FILL_X|PACK_UNIFORM_WIDTH);
  
  new FXButton(addRemovFrame,"&More\tAdd a new search condition", NULL, this, ID_MORE);
  new FXButton(addRemovFrame,"F&ewer\tRemove last search condition", NULL, this, ID_FEWER);
  
  new FXHorizontalSeparator(this);
  
  int  numOper= 6;
  FXString relOper[] = {"<",">","=","<>","<=",">="};
  m_operators.assign(relOper, &relOper[numOper-1]);
  
  m_factory = new ColWidgetFactory();
  
}


long QueryFrame::onCmdMore(FXObject*,FXSelector,void*)
{
  FXComboBox *temp;
  temp = new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_STATIC|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN); 
  temp->appendItem("AND");
  temp->appendItem("OR");
  temp->setNumVisible(2);
  temp->create();
  
  temp = new FXComboBox(m_searchFrame, 0, this, ID_COLSELECT, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_STATIC|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
  int i;
  FXComboBox *firstColSel = (FXComboBox *) m_searchFrame->getFirst();
  for (i = 0; i < firstColSel->getNumItems(); i++)
    temp->appendItem(firstColSel->getItemText(i), firstColSel->getItemData(i));
  temp->setNumVisible(10);
  temp->create();
  
  temp = new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_STATIC|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);  
  for (i = 0; i < m_operators.size(); i++)
    temp->appendItem(m_operators[i]);
  temp->setNumVisible(6);
  temp->create();
    
  FXHorizontalFrame *tempFrame = new FXHorizontalFrame(m_searchFrame,LAYOUT_FILL_X|LAYOUT_FILL_Y|LAYOUT_FILL_COLUMN,  
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  tempFrame->setBackColor(getApp()->getBackColor());
  tempFrame->create();
  temp = new FXComboBox(tempFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_NORMAL|LAYOUT_FILL_X); 
  temp->create();
  m_searchFrame->recalc();
  return 1;
}


long QueryFrame::onCmdFewer(FXObject*,FXSelector,void*)
{
  FXWindow *temp;
  if (m_searchFrame->getNumRows() == 1)
    return 1;
  for (int j = 0; j < m_searchFrame->getNumColumns()-1; j++)
    {
      temp = m_searchFrame->getLast();
      temp->destroy();
      delete temp;
    }
  temp = m_searchFrame->getLast();
  temp->destroy();
  delete temp;
  m_searchFrame->recalc();
  return 1;
}


void QueryFrame::updateColumnSelection(const FXCheckList *colList)
{
  FXWindow *temp; 
  while (m_searchFrame->getNumRows() > 1)
    {
      onCmdFewer(NULL,0,NULL);
    }
    
  FXComboBox *colSelect = (FXComboBox *) m_searchFrame->getFirst();
  colSelect->clearItems();
  
  int i;
  for (i = 0; i < colList->getNumItems(); i++)
    colSelect->appendItem(colList->getItemText(i), colList->getItemData(i));
  colSelect->setNumVisible(10);

  FXComboBox *operSelect = (FXComboBox *) m_searchFrame->childAtRowCol(0,1);
  for (i = 0; i < m_operators.size(); i++)
    operSelect->appendItem(m_operators[i]);
  operSelect->setNumVisible(6);
  
  if (rdbModel::Column* col = (rdbModel::Column *)colList->getItemData(colList->getCurrentItem()))
    {
      FXHorizontalFrame *colWidgetFrame = (FXHorizontalFrame *)m_searchFrame->childAtRowCol(0,2);
      FXWindow *child = colWidgetFrame->getFirst();
      child->destroy();
      delete child;
      m_factory->createColWidget(colWidgetFrame, col);
      m_searchFrame->recalc();
    }
}

long QueryFrame::onSelectCol(FXObject *sender, FXSelector, void*)
{
  FXComboBox *colList = (FXComboBox *) sender;
  FXint row, col;
  row = m_searchFrame->rowOfChild(colList);
  
  if (rdbModel::Column* col = (rdbModel::Column *)colList->getItemData(colList->getCurrentItem()))
  {
    FXHorizontalFrame *colWidgetFrame = (FXHorizontalFrame *)m_searchFrame->childAtRowCol(row,2);
    FXWindow *child = colWidgetFrame->getFirst();
    child->destroy();
    delete child;
    m_factory->createColWidget(colWidgetFrame, col);
    m_searchFrame->recalc();
  }
  return 1;
}


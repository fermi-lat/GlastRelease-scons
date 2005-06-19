
#include "QueryFrame.h"
#include "FXCheckList.h"
#include "ColWidgetFactory.h"
#include "RdbGuiWindow.h"

#include "rdbModel/Tables/Assertion.h"

#include <string>

// Message Map TableColumnList class
FXDEFMAP(QueryFrame) QueryFrameMap[]={

  //__Message_Type_____________ID________________________Message_Handler_____
  FXMAPFUNC(SEL_COMMAND,   QueryFrame::ID_MORE,          QueryFrame::onCmdMore),
  FXMAPFUNC(SEL_COMMAND,   QueryFrame::ID_FEWER,         QueryFrame::onCmdFewer),
  FXMAPFUNC(SEL_COMMAND,   QueryFrame::ID_COLSELECT,     QueryFrame::onSelectCol),
  FXMAPFUNC(SEL_COMMAND,   QueryFrame::ID_QUERY,         QueryFrame::onQuery)
  };

// Object implementation
FXIMPLEMENT(QueryFrame,FXVerticalFrame,QueryFrameMap,ARRAYNUMBER(QueryFrameMap))



QueryFrame::QueryFrame(FXComposite *owner, RdbGUIWindow *target):
  FXVerticalFrame(owner, LAYOUT_FILL_X|LAYOUT_FILL_Y, 0, 0, 0, 0, 0, 0, 0, 0),
  m_target(target)
{

  // Search frame scroll window
  FXScrollWindow *srcScrollWindow = new FXScrollWindow(this, LAYOUT_FILL_X|LAYOUT_FILL_Y);
  
  // Matrix of FXComboBox 
  m_searchFrame = new FXMatrix(srcScrollWindow, 4, MATRIX_BY_COLUMNS|LAYOUT_FILL_X|
      LAYOUT_FILL_Y|FRAME_SUNKEN|FRAME_THICK|PACK_UNIFORM_WIDTH);
  m_searchFrame->setBackColor(getApp()->getBackColor());
  
  new FXComboBox(m_searchFrame, 0, this, ID_COLSELECT, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_STATIC|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
  new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_STATIC|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
  FXHorizontalFrame *tempFrame = new FXHorizontalFrame(m_searchFrame,LAYOUT_FILL_X|LAYOUT_FILL_Y|LAYOUT_FILL_COLUMN,  
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  tempFrame->setBackColor(getApp()->getBackColor());
  new FXComboBox(tempFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_NORMAL|LAYOUT_FILL_X);
  
  FXHorizontalFrame *addRemovFrame = new FXHorizontalFrame(this, LAYOUT_FILL_X|PACK_UNIFORM_WIDTH);
  
  m_moreBtn = new FXButton(addRemovFrame,"&More\tAdd a new search condition", NULL, this, ID_MORE);
  m_fewerBtn = new FXButton(addRemovFrame,"F&ewer\tRemove last search condition", NULL, this, ID_FEWER);
  m_sendBtn = new FXButton(addRemovFrame,"&Send\tSend query to the database", NULL, this, 
      ID_QUERY, BUTTON_NORMAL|LAYOUT_RIGHT);
      
  m_moreBtn->disable();
  
  new FXHorizontalSeparator(this);
  
  int  numOper= 6;
  FXString relOper[] = {">", "<", "=","<>","<=",">="};
  m_operators.assign(relOper, &relOper[numOper]);
  
  m_factory = new ColWidgetFactory();
  
  m_tableName = "";
  m_connect = NULL;
  m_queryResult = 0;
  
}


long QueryFrame::onCmdMore(FXObject*,FXSelector,void*)
{
  FXComboBox *temp;
  temp = new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_STATIC|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN); 
  temp->appendItem("AND");
  temp->appendItem("OR");
  temp->setNumVisible(2);
  temp->create();
  
  FXComboBox *colSelect = new FXComboBox(m_searchFrame, 0, this, ID_COLSELECT, 
      FRAME_SUNKEN|FRAME_THICK|COMBOBOX_STATIC|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
  unsigned int i;
  FXComboBox *firstColSel = (FXComboBox *) m_searchFrame->getFirst();
  for (i = 0; i < (unsigned int)firstColSel->getNumItems(); i++)
    colSelect->appendItem(firstColSel->getItemText(i), firstColSel->getItemData(i));
  colSelect->setNumVisible(10);
  colSelect->create();
  
  temp = new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_STATIC|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);  
  for (i = 0; i < m_operators.size(); i++)
    temp->appendItem(m_operators[i]);
  temp->setNumVisible(6);
  temp->create();
    
  FXHorizontalFrame *tempFrame = new FXHorizontalFrame(m_searchFrame,LAYOUT_FILL_X|LAYOUT_FILL_Y|LAYOUT_FILL_COLUMN,  
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  tempFrame->setBackColor(getApp()->getBackColor());
  tempFrame->create();
  if (rdbModel::Column* col = (rdbModel::Column *)colSelect->getItemData(colSelect->getCurrentItem()))
    {
       m_widgets.push_back(m_factory->createColWidget(tempFrame, col));
    }
  else
    {
      temp = new FXComboBox(tempFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_NORMAL|LAYOUT_FILL_X); 
      temp->create();
    }
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
  delete m_widgets[m_widgets.size()-1];
  m_widgets.pop_back();
  
  temp = m_searchFrame->getLast();
  temp->destroy();
  delete temp;
  m_searchFrame->recalc();
  return 1;
}


void QueryFrame::updateColumnSelection(const FXList *tableList, const FXCheckList *colList)
{
  unsigned int i;
  
  m_tableName = tableList->getItemText(tableList->getCurrentItem()).text();
  
  while (m_searchFrame->getNumRows() > 1)
    {
      onCmdFewer(NULL,0,NULL);
    }
    
  for (i = 0; i < m_widgets.size(); i++)
    {
      delete m_widgets[i];
    }
  m_widgets.clear();
    
  FXComboBox *colSelect = (FXComboBox *) m_searchFrame->getFirst();
  colSelect->clearItems();
  
  for (i = 0; i < (unsigned int) colList->getNumItems(); i++)
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
      m_widgets.push_back(m_factory->createColWidget(colWidgetFrame, col));
      m_searchFrame->recalc();
    }
}

long QueryFrame::onSelectCol(FXObject *sender, FXSelector, void*)
{
  FXComboBox *colList = (FXComboBox *) sender;
  FXint row;
  row = m_searchFrame->rowOfChild(colList);
  
  if (rdbModel::Column* col = (rdbModel::Column *)colList->getItemData(colList->getCurrentItem()))
  {
    FXHorizontalFrame *colWidgetFrame = (FXHorizontalFrame *)m_searchFrame->childAtRowCol(row,2);
    FXWindow *child = colWidgetFrame->getFirst();
    child->destroy();
    delete child;
    delete m_widgets[row];
    m_widgets[row] = m_factory->createColWidget(colWidgetFrame, col);
    m_searchFrame->recalc();
  }
  return 1;
}


long QueryFrame::onQuery(FXObject*,FXSelector,void*)
{

  rdbModel::Assertion::Operator *whereOp = buildOperator(0);
  
  rdbModel::Assertion *where= new rdbModel::Assertion(whereOp);
  
  if (m_connect)
    {
      int i;
      std::vector<std::string> getCols;
      std::vector<std::string> orderCols;   // needed  by the select interace (I just leave it empty)
      
      FXComboBox *firstColSel = (FXComboBox *) m_searchFrame->getFirst();
      for (i = 0; i < firstColSel->getNumItems(); i++)
      {
        getCols.push_back(firstColSel->getItemText(i).text());
      }

      if (m_target->getPrimaryName() != "")
        orderCols.push_back(m_target->getPrimaryName());
       
      if (m_queryResult) 
        {
          delete m_queryResult;
          m_queryResult = NULL;
        }
      m_queryResult = m_connect->select(m_tableName, getCols, orderCols, where);
      return m_target && m_target->handle(NULL,FXSEL(SEL_COMMAND,ID_QUERY),NULL);
    }  

  return 0;
}


rdbModel::Assertion::Operator* QueryFrame::buildCompOperator(std::string col, 
    std::string comp, std::string value)
{
  rdbModel::OPTYPE opType;
  
  if (comp == ">") 
    opType = rdbModel::OPTYPEgreaterThan;  
  else if (comp == "<")
    opType = rdbModel::OPTYPElessThan;
  else if  (comp ==  "=")
    opType = rdbModel::OPTYPEequal;
  else if (comp == "<>")
    opType = rdbModel::OPTYPEnotEqual;
  else if (comp == "<=")
    opType = rdbModel::OPTYPElessOrEqual;
  else if (comp == ">=")
    opType = rdbModel::OPTYPEgreaterOrEqual;
    
  return new rdbModel::Assertion::Operator(opType, col, value, 
                                           rdbModel::FIELDTYPEold,
                                           rdbModel::FIELDTYPElit);
}

rdbModel::Assertion::Operator* QueryFrame::buildOperator(int row)
{
  rdbModel::Assertion::Operator *leftOper;
  std::vector<rdbModel::Assertion::Operator* > children;
  
  FXComboBox *column, *compOp, *conjOp;
  std::string leftArg, comp, rightArg;
  

  column = (FXComboBox *) m_searchFrame->childAtRowCol(row,0);
  compOp = (FXComboBox *) m_searchFrame->childAtRowCol(row,1);
  
  leftArg = column->getText().text();
  comp = compOp->getText().text();
  rightArg = m_widgets[row]->getValue();
  
  children.push_back(buildCompOperator(leftArg, comp, rightArg));
  
  while (row + 1 < m_searchFrame->getNumRows())
    {
      conjOp = (FXComboBox *)m_searchFrame->childAtRowCol(row,3);

      if (conjOp->getText()=="AND")
        {
          row++;
          column = (FXComboBox *) m_searchFrame->childAtRowCol(row,0);
          compOp = (FXComboBox *) m_searchFrame->childAtRowCol(row,1);
          
          leftArg = column->getText().text();
          comp = compOp->getText().text();
          rightArg = m_widgets[row]->getValue();
          children.push_back(buildCompOperator(leftArg, comp, rightArg));
        }
      else
        break;
    }
    
  if (children.size() > 1)
    {
      leftOper = new rdbModel::Assertion::Operator(rdbModel::OPTYPEand, children);
    }
  else
    {
      leftOper = children[0];
    }
    
  children.clear();
  if (row + 1 < m_searchFrame->getNumRows())
    {
      children.push_back(leftOper);
      children.push_back(buildOperator(row + 1));
      leftOper = new rdbModel::Assertion::Operator(rdbModel::OPTYPEor, children);
    }
    
  return leftOper;
}


void QueryFrame::reset()
{
  m_connect = NULL;
}

void QueryFrame::setEnabled(bool flag)
{
  if (flag)
    {
      m_moreBtn->enable();
      m_fewerBtn->enable();
      m_sendBtn->enable();
    }
  else
    {
      m_moreBtn->disable();
      m_fewerBtn->disable();
      m_sendBtn->disable();
    }
}

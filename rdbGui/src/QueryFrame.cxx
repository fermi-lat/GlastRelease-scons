
#include "QueryFrame.h"

// Message Map TableColumnList class
FXDEFMAP(QueryFrame) QueryFrameMap[]={

  //__Message_Type_____________ID________________________Message_Handler_____
  FXMAPFUNC(SEL_COMMAND,   QueryFrame::ID_MORE,          QueryFrame::onCmdMore),
  FXMAPFUNC(SEL_COMMAND,   QueryFrame::ID_FEWER,          QueryFrame::onCmdFewer)
  
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
      LAYOUT_FILL_Y|FRAME_SUNKEN|FRAME_THICK);
  m_searchFrame->setBackColor(FXRGB(255,255,255));
  
  new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_NORMAL|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
  new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_NORMAL|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
  new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_NORMAL|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
  (new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_NORMAL|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN))->hide();
  
  FXHorizontalFrame *addRemovFrame = new FXHorizontalFrame(this, LAYOUT_FILL_X|PACK_UNIFORM_WIDTH);
  
  new FXButton(addRemovFrame,"&More\tAdd a new search condition", NULL, this, ID_MORE);
  new FXButton(addRemovFrame,"F&ewer\tRemove last search condition", NULL, this, ID_FEWER);
  
  new FXHorizontalSeparator(this);
}


long QueryFrame::onCmdMore(FXObject*,FXSelector,void*)
{
  (m_searchFrame->childAtRowCol(m_searchFrame->getNumRows()-1, m_searchFrame->getNumColumns()-1))->show();
  FXComboBox *temp;
  temp = new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_NORMAL|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);
  temp->create();
  temp = new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_NORMAL|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN);  
  temp->create();
  temp = new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_NORMAL|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN); 
  temp->create();
  temp = new FXComboBox(m_searchFrame, 0, NULL, 0, FRAME_SUNKEN|FRAME_THICK|COMBOBOX_NORMAL|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN); 
  temp->create();
  temp->hide();
  m_searchFrame->recalc();
  return 1;
}


long QueryFrame::onCmdFewer(FXObject*,FXSelector,void*)
{
  if (m_searchFrame->getNumRows() == 1)
    return 1;
  for (int j = 0; j < m_searchFrame->getNumColumns(); j++)
    {
      FXWindow *temp = m_searchFrame->getLast();
      temp->destroy();
      delete temp;
    }
  m_searchFrame->getLast()->hide();
  m_searchFrame->recalc();
  return 1;
}
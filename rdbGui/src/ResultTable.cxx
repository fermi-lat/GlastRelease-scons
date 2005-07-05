
#include "ResultTable.h"

#include "fxkeys.h"

#include <vector>

// Message Map ResultTable class
FXDEFMAP(ResultTable) ResultTableMap[]={

  //__Message_Type___________________________ID_______________________________________Message_Handler_____
  FXMAPFUNC(SEL_KEYPRESS,            0,                                    ResultTable::onKeyPress),
  FXMAPFUNC(SEL_COMMAND,             FXTable::ID_COPY_SEL,                 ResultTable::onCmdCopySel),
  FXMAPFUNC(SEL_RIGHTBUTTONRELEASE,  0,                                    ResultTable::onCmdMenuPane),
  FXMAPFUNC(SEL_COMMAND,             ResultTable::ID_UPDATEROW,            ResultTable::onUpdRow),
  FXMAPFUNC(SEL_COMMAND,             ResultTable::ID_COPYROW,              ResultTable::onCopyRow)    
};



// Macro for the GLViewWindow class hierarchy implementation
FXIMPLEMENT(ResultTable,FXTable,ResultTableMap,ARRAYNUMBER(ResultTableMap))


/*******************************************************************************/

ResultTable::ResultTable(FXComposite *p, FXObject* tgt,
                         FXSelector sel, FXuint opts,FXint x,FXint y,FXint w,
                         FXint h,FXint pl,FXint pr,FXint pt,FXint pb):
  FXTable(p, tgt, sel, opts, x, y, w, h, pl, pr, pt, pb)
{
  m_recordActions = new FXMenuPane(this);
  m_updRow = new FXMenuCommand(m_recordActions,"Update Row",NULL, this, ResultTable::ID_UPDATEROW);
  m_copyRow = new FXMenuCommand(m_recordActions,"Copy Row",NULL, this, ResultTable::ID_COPYROW);

}

 
// Set table background colors
void ResultTable::setBackColors(FXColor c00, FXColor c01,FXColor c10, FXColor c11)
{
  setCellColor(0, 0, c00);
  setCellColor(0, 1, c01);
  setCellColor(1, 0, c10);
  setCellColor(1, 1, c11);
}

void ResultTable::showColumn(FXint index)
{
  enableColumn(index);
  formatCol(index);
}

void ResultTable::hideColumn(FXint index)
{
  if ((index >= 0) && (index < ncols) )
    {
      setColumnWidth(index, 0);
      disableColumn(index);
    }
}

void ResultTable::disableColumn(FXint index)
{
  int r;
  if ((index >= 0) && (index < ncols) )
    {
      for (r = 0; r < nrows; r++)
        {
          getItem(r, index)->setEnabled(false);
        }
    }  
}

void ResultTable::enableColumn(FXint index)
{
  int r;
  if ((index >= 0) && (index < ncols) )
    {
      for (r = 0; r < nrows; r++)
        {
          getItem(r, index)->setEnabled(true);
        }
    }  
}

void ResultTable::formatCol(FXint index)
{
  if ((index >= 0) && (index < ncols))
    {
      int i;
      int colWidth = getDefColumnWidth();
      for (i = 0; i < nrows; i++)
        {
          setItemJustify(i,index, FXTableItem::LEFT);
          int currWidth = getItem(i,index)->getWidth(this);
          colWidth = (colWidth < currWidth) ? currWidth : colWidth;
        }
        
      setColumnWidth(index, colWidth);
    }
}


void ResultTable::format()
{
  int i,j;

  // for each column, find the maximum width
  
  std::vector<int> colWidth(ncols,0);
  
  int rowWidth = 0;
  for (i = 0; i < nrows; i++)
    {
      int currWidth = getRowText(i).length()*getFont()->getFontWidth();
      rowWidth = (rowWidth < currWidth) ? currWidth : rowWidth;
    }
  setRowHeaderWidth(rowWidth);
      

  for (j = 0; j < ncols; j++)
    {
      formatCol(j);
    }

}


long ResultTable::onKeyPress(FXObject*,FXSelector,void* ptr)
{
  FXEvent *event = (FXEvent*) ptr;
  if(event->state&CONTROLMASK)
    {
      if (event->code == KEY_c)
        {
          handle(this,FXSEL(SEL_COMMAND,FXTable::ID_COPY_SEL),NULL);
          return 1;
        }
    }
  return 0;
}  


// Extract cells from given range as text
void ResultTable::extractText(FXchar*& text,FXint& size,FXint startrow,FXint endrow,FXint startcol,FXint endcol,FXString cs,FXchar rs) const 
{
  register FXchar *ptr;
  register FXuint sz=0;
  register FXint r,c;
  FXString string;

  // Verify range
  if(startrow<0 || startcol<0 || nrows<=endrow || ncols<=endcol){ fxerror("%s::extractText: index out of range.\n",getClassName()); }

  // Initialize
  text=NULL;
  size=0;
  int colSepSize = cs.length();

  // Non-empty range
  if(startrow<=endrow && startcol<=endcol){
    for(r=startrow; r<=endrow; r++){
      for(c=startcol; c<=endcol; c++){
        if (getItem(r,c)->isEnabled())
          sz+=getItemText(r,c).length() + colSepSize;
        }
      }
    if(FXMALLOC(&text,FXchar,sz+1)){
      size=sz;
      ptr=text;
      for(r=startrow; r<=endrow; r++){
        for(c=startcol; c<=endcol; c++){
          if (getItem(r,c)->isEnabled())
            {
              string=getItemText(r,c);
              memcpy(ptr,string.text(),string.length());
              ptr+=string.length();
              if (c == endcol)
                {
                  *ptr++ = rs;
                }
              else
                {
                  memcpy(ptr, cs.text(), colSepSize);
                  ptr+=colSepSize;
                }
            }
          }
        }
      *ptr='\0';        // Its there but not accounted for...
      }
    }
}

// Copy selection
long ResultTable::onCmdCopySel(FXObject*,FXSelector,void*){
  if(isAnythingSelected()){
    FXDragType types[3];
    types[0]=stringType;
    types[1]=textType;
    types[2]=csvType;
    if(acquireClipboard(types,3)){
      FXFREE(&clipbuffer);
      extractText(clipbuffer,cliplength,selection.fm.row,selection.to.row,selection.fm.col,selection.to.col);
      }
    }
  return 1;
  }
  
long ResultTable::onCmdMenuPane(FXObject* sender,FXSelector sel,void* ptr)
{
  FXTable::onRightBtnRelease(sender,sel,ptr);
  FXEvent *event = (FXEvent*) ptr;
  m_recordActions->popup(NULL,event->root_x,event->root_y);
  FXint x,y;
  FXuint buttons;
  getCursorPosition(x,y,buttons);
  m_selRow = rowAtY(y);
  if ((m_selRow>=0) && (m_selRow < nrows) && (getItemText(0, 0) !="No data") &&
      (getItemText(0, 0) != "Null query result") )
    {
      selectRow(m_selRow);
      m_updRow->enable();
    }
  else
    {
      m_updRow->disable();
    }
  return 1;
}

long ResultTable::onCopyRow(FXObject*,FXSelector,void*)
{
  return target && target->handle(NULL,FXSEL(SEL_COMMAND,ID_COPYROW),(void *)m_selRow);
}

long ResultTable::onUpdRow(FXObject*,FXSelector,void*)
{
  return target && target->handle(NULL,FXSEL(SEL_COMMAND,ID_UPDATEROW),(void *)m_selRow);
}

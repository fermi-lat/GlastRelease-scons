
#include "ResultTable.h"

#include <vector>

ResultTable::ResultTable(FXComposite *p, FXObject* tgt,
                         FXSelector sel, FXuint opts,FXint x,FXint y,FXint w,
                         FXint h,FXint pl,FXint pr,FXint pt,FXint pb):
  FXTable(p, tgt, sel, opts, x, y, w, h, pl, pr, pt, pb)
{

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
  formatCol(index);
}

void ResultTable::hideColumn(FXint index)
{
  if ((index >= 0) && (index < ncols) )
    setColumnWidth(index, 0);
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

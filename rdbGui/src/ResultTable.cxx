
#include "ResultTable.h"

#include <vector>


ResultTable::ResultTable(FXComposite *p, FXObject* tgt,
                         FXSelector sel, FXuint opts,FXint x,FXint y,FXint w,
                         FXint h,FXint pl,FXint pr,FXint pt,FXint pb):
  FXTable(p, tgt, sel, opts, x, y, w, h, pl, pr, pt, pb)
{

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
    colWidth[j] = getColumnWidth(j);
  
  for (i = 0; i < nrows; i++)
    for (j = 0; j < ncols; j++)
      {
        setItemJustify(i,j, FXTableItem::LEFT);
        int currWidth = getItem(i,j)->getWidth(this);
        colWidth[j] = (colWidth[j] < currWidth) ? currWidth : colWidth[j];
      }

  for (j = 0; j < ncols; j++)
    setColumnWidth(j, colWidth[j]);
}  


#include "ResultTable.h"


ResultTable::ResultTable(FXComposite *p, FXObject* tgt,
                         FXSelector sel, FXuint opts,FXint x,FXint y,FXint w,
                         FXint h,FXint pl,FXint pr,FXint pt,FXint pb):
  FXTable(p, tgt, sel, opts, x, y, w, h, pl, pr, pt, pb)
{

}


void ResultTable::format()
{
  int i,j;

  // Return if table is empty, otherwise store column and row boundary
  if (ncols < 1 || nrows < 1)
    return;
  
  //  setLeadingColumns(1);
  //  setLeadingRows(1);
  
  for(j = 0; j < ncols; j++)
    {
      //getItem(0,i)->setButton(true);
      getItem(0,i)->setJustify(FXTableItem::TOP); 
    }
    
//   for(i = 0; i < nrows; i++)
//     {
//       getItem(i,0)->setButton(true);
//     }
}  

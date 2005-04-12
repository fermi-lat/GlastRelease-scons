
#ifndef RESULTTABLE_H
#define RESULTTABLE_H

#include "fx.h"
#include <string>

class ResultTable: public FXTable
{
  FXDECLARE(ResultTable) 
   
  public:
  
    enum {
      ID_UPDATEROW = FXTable::ID_LAST,
      ID_COPYROW
    };
  
    ResultTable(FXComposite *p, FXObject* tgt=NULL,FXSelector sel=0,
                FXuint opts=0,FXint x=0,FXint y=0,FXint w=0,FXint h=0,
                FXint pl=DEFAULT_MARGIN,FXint pr=DEFAULT_MARGIN,
                FXint pt=DEFAULT_MARGIN,FXint pb=DEFAULT_MARGIN);
       
    void setBackColors(FXColor, FXColor, FXColor, FXColor);   
    void formatCol(FXint);     
    void format();
    void showColumn(FXint);
    void hideColumn(FXint);  
    void enableColumn(FXint);
    void disableColumn(FXint);
    
    void setTableName(std::string tbName){m_tableName = tbName;}
    std::string getTableName(){return m_tableName;}
    
    /// Extract cells from given range as text.
    void extractText(FXchar*& text,FXint& size,FXint startrow,FXint endrow,FXint startcol,FXint endcol,FXString cs="\t",FXchar rs='\n') const;
    
    long onKeyPress(FXObject*,FXSelector,void*);
    long onCmdCopySel(FXObject*,FXSelector,void*);
    long onCmdMenuPane(FXObject*,FXSelector,void*);
    long onUpdRow(FXObject*,FXSelector,void*);
    long onCopyRow(FXObject*,FXSelector,void*);
    
  protected:
    ResultTable() {};
    ResultTable(const ResultTable&);
    
  private:
    std::string m_tableName;
    FXMenuPane *m_recordActions;
    FXMenuCommand *m_updRow;
    FXMenuCommand *m_copyRow;
    FXint m_selRow;
    
};


#endif // end of RESULTTABLE_H

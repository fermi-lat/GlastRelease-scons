
#ifndef RESULTTABLE_H
#define RESULTTABLE_H

#include "fx.h"

class ResultTable: public FXTable
{
  public:
  
    ResultTable(FXComposite *p, FXObject* tgt=NULL,FXSelector sel=0,
                FXuint opts=0,FXint x=0,FXint y=0,FXint w=0,FXint h=0,
                FXint pl=DEFAULT_MARGIN,FXint pr=DEFAULT_MARGIN,
                FXint pt=DEFAULT_MARGIN,FXint pb=DEFAULT_MARGIN);
                
    void format();  
    
  protected:
    ResultTable() {};
    ResultTable(ResultTable&) {};
};


#endif // end of RESULTTABLE_H


#ifndef QUERYFRAME_H
#define QUERYFRAME_H


#include "fx.h"

class QueryFrame: public FXVerticalFrame
{
  FXDECLARE(QueryFrame)
  
 public:
 
  enum{
    ID_QUERYFRAME=FXVerticalFrame::ID_LAST,
    ID_MORE,
    ID_FEWER
  }; 
  
  QueryFrame(FXComposite *);
  
  long onCmdMore(FXObject*,FXSelector,void*);
  long onCmdFewer(FXObject*,FXSelector,void*);
  
 protected:
  QueryFrame(){}
  QueryFrame(QueryFrame&){} 
  
 private:
  FXMatrix *m_searchFrame;                  // Martix of FXComboBox containing search conditions


};

#endif
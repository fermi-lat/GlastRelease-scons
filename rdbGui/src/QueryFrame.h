
#ifndef QUERYFRAME_H
#define QUERYFRAME_H


#include "fx.h"

#include <vector>
#include "FXCheckList.h"

class ColWidgetFactory;

class QueryFrame: public FXVerticalFrame
{
  FXDECLARE(QueryFrame)
  
 public:
 
  enum{
    ID_QUERYFRAME=FXVerticalFrame::ID_LAST,
    ID_MORE,
    ID_FEWER,
    ID_COLSELECT
  }; 
  
  QueryFrame(FXComposite *);
  
  long onCmdMore(FXObject*,FXSelector,void*);
  long onCmdFewer(FXObject*,FXSelector,void*);
  long onSelectCol(FXObject *sender, FXSelector, void*);
  
  void updateColumnSelection(const FXCheckList *colList);
  
 protected:
  QueryFrame(){}
  QueryFrame(QueryFrame&){} 
  
 private:
  FXMatrix *m_searchFrame;                  // Martix of FXComboBox containing search conditions
  std::vector<FXString> m_operators;     // vector of comparison operators
  ColWidgetFactory* m_factory;  


};

#endif
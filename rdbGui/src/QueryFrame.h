
#ifndef QUERYFRAME_H
#define QUERYFRAME_H


#include "fx.h"

#include <vector>
#include "FXCheckList.h"

#include "rdbModel/Db/Connection.h"
#include "rdbModel/Tables/Assertion.h"
#include "rdbModel/Db/ResultHandle.h"

class ColWidget;
class ColWidgetFactory;
class RdbGUIWindow;

class QueryFrame: public FXVerticalFrame
{
  FXDECLARE(QueryFrame)
  
 public:
 
  enum{
    ID_QUERYFRAME=FXVerticalFrame::ID_LAST,
    ID_MORE,
    ID_FEWER,
    ID_COLSELECT,
    ID_QUERY
  }; 
  
  QueryFrame(FXComposite *, RdbGUIWindow *target = NULL);
  
  long onCmdMore(FXObject*,FXSelector,void*);
  long onCmdFewer(FXObject*,FXSelector,void*);
  long onSelectCol(FXObject *sender, FXSelector, void*);
  long onQuery(FXObject *sender, FXSelector, void*);
  
  void updateColumnSelection(const FXList *tableList, const FXCheckList *colList);
  void setConnection(rdbModel::Connection* con){m_connect = con;}
  rdbModel::ResultHandle* getQueryResult() const {return m_queryResult;}
  void reset();
  void setEnabled(bool);

  
 protected:
  QueryFrame(){}
  QueryFrame(const QueryFrame&);
  
 private:
  RdbGUIWindow *m_target;                    // The target of some messages sent by this widget
  FXMatrix *m_searchFrame;               // Martix of FXComboBox containing search conditions
  std::vector<FXString> m_operators;     // vector of comparison operators
  ColWidgetFactory* m_factory;  
  std::vector<ColWidget*> m_widgets;
  rdbModel::Connection* m_connect;       // pointer to the DB connection class
  std::string m_tableName;               // Name of the selected table
  rdbModel::ResultHandle *m_queryResult; // object containing the result of a query
  FXButton* m_moreBtn;                   // More button  
  FXButton* m_fewerBtn;                  // Fewer button
  FXButton* m_sendBtn;                   // Send button

  rdbModel::Assertion::Operator* QueryFrame::buildCompOperator(std::string col, 
      std::string comp, std::string value);
  rdbModel::Assertion::Operator* buildOperator(int row);
  
};

#endif

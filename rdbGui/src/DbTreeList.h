
#ifndef DBTREELIST_H
#define DBTREELIST_H

#include "fx.h"

#include "rdbModel/Management/Visitor.h"
#include "rdbModel/Rdb.h"
#include "rdbModel/Tables/Table.h"
#include "rdbModel/Tables/Assertion.h"
#include "rdbModel/Tables/Column.h"

#include <vector>

namespace rdbModel {
  class InsertRow;
  class Supersede;
  class InterRow;
  class Query;
  class Set;
}

class DbTreeList : public FXTreeList, public rdbModel::Visitor 
{

 private: 
  FXint m_lastx, m_lasty;
  std::vector<FXTreeItem*> m_rootAtLevel;  
  FXICOIcon *dbIcon, *tbIcon;

 public:
  
  DbTreeList(FXComposite *p, FXObject* tgt=NULL,FXSelector sel=0,
	     FXuint opts=TREELIST_NORMAL,FXint x=0,FXint y=0,FXint w=0,FXint h=0);
             
  void init();

  FXint last_x() {return m_lastx;}
  FXint last_y() {return m_lasty;}

  void rememberPos(FXint x, FXint y);


  /// Recursively scan the whole tree and store all items into the vector
  void scanForAll(const FXTreeItem *subTree, std::vector<FXString> *all);

  /// Recursively scan the whole tree and store selected items into the array  
  void scanForSelected(const FXTreeItem *subTree, std::vector<FXString> *sel);

  /// Get the array of selected items (wrapper for scanForSelected)
  std::vector<FXString>* getSelectedText();
  
  /// Remove subtree possibly including the subtree root
  void clearSubTree(FXTreeItem *root, FXbool inclusive = false);
  
 
 private:
  
  /// Visitors to build the tree widget representing the database     
  rdbModel::Visitor::VisitorState visitRdb(rdbModel::Rdb *rdb);
  rdbModel::Visitor::VisitorState visitTable(rdbModel::Table *table);
  rdbModel::Visitor::VisitorState visitColumn(rdbModel::Column *column);
  rdbModel::Visitor::VisitorState visitIndex(rdbModel::Index *index);
  rdbModel::Visitor::VisitorState visitAssertion(rdbModel::Assertion *assertion);

  rdbModel::Visitor::VisitorState visitInsertNew(rdbModel::InsertNew*) {
    return Visitor::VCONTINUE;
  }
  rdbModel::Visitor::VisitorState visitSupersede(rdbModel::Supersede* ) {
    return Visitor::VCONTINUE;
  }
  rdbModel::Visitor::VisitorState visitQuery(rdbModel::Query* ){
    return Visitor::VCONTINUE;
  }
  rdbModel::Visitor::VisitorState visitInterRow(rdbModel::InterRow* ) {
    return Visitor::VCONTINUE;
  }
  rdbModel::Visitor::VisitorState visitSet(rdbModel::Set*) {
    return Visitor::VCONTINUE;
  }


};


#endif // end of DBTREELIST_H

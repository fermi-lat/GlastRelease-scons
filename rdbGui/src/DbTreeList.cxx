
#include "DbTreeList.h"
#include "Icons.h"

#include <string>


DbTreeList::DbTreeList(FXComposite *p, FXObject* tgt,FXSelector sel,
		       FXuint opts,FXint x,FXint y,FXint w,FXint h):
  FXTreeList(p, tgt, sel, opts, x, y, w, h),
  m_rootAtLevel(0) 
{
  dbIcon = new FXICOIcon(getApp(),dbIconImg);
  dbIcon->create();
  tbIcon = new FXICOIcon(getApp(),tbIconImg);
  tbIcon->create();
}

void DbTreeList::init()
{
  clearItems();
  m_rootAtLevel.clear();
}

void DbTreeList::rememberPos(FXint x, FXint y)
{
  m_lastx = x;
  m_lasty = y;
}


// Recursively scan the whole tree and store all items into the array
void DbTreeList::scanForAll(const FXTreeItem *subTree, std::vector<FXString> *all)
{
  const FXTreeItem *item = subTree->getFirst();
  for (int i = 0; i < subTree->getNumChildren(); i++)
    {
      all->push_back(getItemText(item));
      scanForAll(item, all);
      item = subTree->getNext();
    }
  
}


// Recursively scan the whole tree and store selected items into the array
void DbTreeList::scanForSelected(const FXTreeItem *subTree, std::vector<FXString> *sel)
{
  const FXTreeItem *item = subTree->getFirst();
  for (int i = 0; i < subTree->getNumChildren(); i++)
    {
      if (isItemSelected(item))
	sel->push_back(getItemText(item));
      if (isItemExpanded(item))
	scanForSelected(item, sel);
      item = subTree->getNext();
    }
    
}


// Get the array of selected items (wrapper for scanForSelected)
std::vector<FXString>* DbTreeList::getSelectedText()
{
  std::vector<FXString> *sel = new std::vector<FXString>();
  scanForSelected(getFirstItem(), sel);
  return sel;
}
  

// Remove subtree possibly including the subtree root
void DbTreeList::clearSubTree(FXTreeItem *root, FXbool inclusive)
{
  removeItems(root->getFirst(), root->getLast());
  if (inclusive)
    removeItem(root);
}


rdbModel::Visitor::VisitorState DbTreeList::visitRdb(rdbModel::Rdb *rdb)
{
  FXTreeItem *root = addItemLast(NULL, rdb->getDbName().c_str(), dbIcon, dbIcon);
  m_rootAtLevel.push_back(root);
  expandTree(root);
  return rdbModel::Visitor::VCONTINUE;
}

rdbModel::Visitor::VisitorState DbTreeList::visitTable(rdbModel::Table *table)
{
  static int level = m_rootAtLevel.size();
  FXTreeItem *item = addItemLast(m_rootAtLevel[level - 1], table->getName().c_str(), tbIcon, tbIcon);
  if (m_rootAtLevel.size() > level)
    m_rootAtLevel[level] = item;
  else
    m_rootAtLevel.push_back(item);
  return rdbModel::Visitor::VCONTINUE;
}

rdbModel::Visitor::VisitorState DbTreeList::visitColumn(rdbModel::Column *column)
{
  static int level = m_rootAtLevel.size();
  FXTreeItem *item = addItemLast(m_rootAtLevel[level - 1], column->getName().c_str());
  if (m_rootAtLevel.size() > level)
    m_rootAtLevel[level] = item;
  else
    m_rootAtLevel.push_back(item);  
  return rdbModel::Visitor::VCONTINUE;
}

rdbModel::Visitor::VisitorState DbTreeList::visitIndex(rdbModel::Index *index)
{
  return rdbModel::Visitor::VCONTINUE;
}

rdbModel::Visitor::VisitorState DbTreeList::visitAssertion(rdbModel::Assertion *assertion)
{
  return rdbModel::Visitor::VCONTINUE;
}

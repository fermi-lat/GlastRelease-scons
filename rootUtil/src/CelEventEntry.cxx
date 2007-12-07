
/*
* Project: GLAST
* Package: rootUtil
*    File: $Id$
* Authors:
*   EC, Eric Charles,    SLAC              echarles@slac.stanford.edu
*
* Copyright (c) 2007
*                   Regents of Stanford University. All rights reserved.
*
*/

#include <rootUtil/CelEventEntry.h>
#include <rootUtil/FileUtil.h>
#include <rootUtil/BgDataHandle.h>
#include <rootUtil/CelFileAndTreeSet.h>

#include <TTree.h>
#include <Riostream.h>



//====================================================================
// Construction and assignment
// One of the use-case is that one could want to manipulate
// two CELs at the same time, on efor reading and one for
// writing. In this context, it could make sense to construct
// first the reading CEL, and duplicate it so to get a
// write one with same structure. This is why we have copy
// constructors which are copyig the BgDataHandle values
// but relinking them to the new BranchGroup. Also, the assignment
// is assigning only the value, not the link between handles and group,
// so that ultimately we could do writeCel = readCel which would
// duplicate values of the current Event.
// As a summary : copy construction is duplicating the structure,
// and assigment only concerns the current event values.
//====================================================================


ClassImp(CelEventEntry) ;

CelEventEntry::CelEventEntry()
 : _componentName(),
   _entryIndex(-1,*this,"Event_EntryIndex"),  
   _treeIndex(FileUtil::NOKEY,*this,"Event_TreeIndex")
 { BgDataHandleInstance::init() ; }
	      
CelEventEntry::CelEventEntry( const TString & componentName )
 : _componentName(componentName),
   _entryIndex(-1,*this,"Event_EntryIndex"),  
   _treeIndex(FileUtil::NOKEY,*this,"Event_TreeIndex")
 { BgDataHandleInstance::init() ; }	      
	      
CelEventEntry::CelEventEntry( const CelEventEntry & other )
 : BranchGroup(),
   _componentName(other._componentName),
   _entryIndex(other._entryIndex,*this,"Event_EntryIndex"),  
   _treeIndex(other._treeIndex,*this,"Event_TreeIndex")
 {}	      

CelEventEntry::~CelEventEntry()
 {}

CelEventEntry & CelEventEntry::operator=( const CelEventEntry & other )
 {
  if ( this == &other ) return *this ;
  _componentName = other.componentName() ; // why not assert it is the same ?
  _entryIndex  = other.entryIndex() ;
  _treeIndex = other.treeIndex() ;
  return *this ;
 }



//====================================================================
// Connection with trees
//====================================================================


// Grab the current status from 'tree'
//
// The event index is set from tree.GetReadEntry();
// The tree index is looked up in handle
void CelEventEntry::set( TTree & tree, CelFileAndTreeSet & handle )
 {
  _entryIndex = tree.GetReadEntry() ;
//  std::cout
//    <<"[CelEventEntry::set] setting entry " <<_entryIndex
//    <<" for tree "<<tree.GetDirectory()->GetUUID().AsString()<<"/"<<tree.GetName()
//    <<" at "<<&tree
//    <<std::endl ;
  UShort_t tIdx = handle.getKey(&tree) ;
  if ( tIdx == FileUtil::NOKEY )
   {
	tIdx = handle.addTree(tree) ;
//    std::cout
//      <<"[CelEventEntry::set] new tree "<<&tree
//      <<" has received index "<<tIdx
//      <<std::endl ;
   }
  _treeIndex = tIdx ;
 }

// Read an event using handle to look up the correct tree
//
// The tree lookup uses _treeIndex to get the tree from 'handle'
// Then the Entry _entryIndex is loaded using LoadTree(_entryIndex);
Int_t CelEventEntry::read( const CelFileAndTreeSet & handle )
 {
  TTree * t = handle.getTree(_treeIndex) ;
  if ( 0 == t ) return -1 ;
  return t->LoadTree(_entryIndex) ;
 }

// Use 'handle' to look up a tree 
//
// The tree lookup uses _treeIndex to get the tree from 'handle'
TTree * CelEventEntry::getTree( const CelFileAndTreeSet & handle ) const
 {
  return handle.getTree(_treeIndex) ;
 }

// Print the information about the current event to cout
// 
// The print format is treeIndex|entryIndex
void CelEventEntry::printInfo() const
 {
  std::cout
    <<setw(2)<<std::right<< _treeIndex
    <<'|'
    <<setw(2)<<std::left<< _entryIndex ;
 }


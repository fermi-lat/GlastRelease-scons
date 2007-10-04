
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

#include "rootUtil/CelEventEntry.h"
#include "rootUtil/FileUtil.h"
#include "rootUtil/BgDataHandle.h"
#include "rootUtil/CelFileAndTreeSet.h"

#include <TTree.h>
#include <Riostream.h>



//====================================================================
// 
//====================================================================


ClassImp(CelEventEntry) ;

CelEventEntry::CelEventEntry()
 : _componentName(),
   _entryIndex(-1,*this,"EntryIndex"),  
   _treeIndex(FileUtil::NOKEY,*this,"TreeIndex")
 { BgDataHandleInstance::init() ; }
	      
CelEventEntry::CelEventEntry( const TString & componentName )
 : _componentName(componentName),
   _entryIndex(-1,*this,"EntryIndex"),  
   _treeIndex(FileUtil::NOKEY,*this,"TreeIndex")
 { BgDataHandleInstance::init() ; }	      
	      
CelEventEntry::CelEventEntry( const CelEventEntry & other)
 : BranchGroup(),
   _componentName(other._componentName),
   _entryIndex(other._entryIndex,*this,"EntryIndex"),  
   _treeIndex(other._treeIndex,*this,"TreeIndex")
 { BgDataHandleInstance::init() ; }	      

// D'tor is a no-op
CelEventEntry::~CelEventEntry()
 {}

// Assignment operator
CelEventEntry & CelEventEntry::operator=( const CelEventEntry & other )
 {
  if ( this == &other ) return *this ;
  _componentName = other.componentName() ;
  _entryIndex  = other.entryIndex() ;
  _treeIndex = other.treeIndex() ;
  return *this ;
 }



//====================================================================
// 
//====================================================================


// Grab the current status from 'tree'
//
// The event index is set from tree.GetReadEntry();
// The tree index is looked up in handle
void CelEventEntry::set( TTree & tree, CelFileAndTreeSet & handle )
 {
  _entryIndex = tree.GetReadEntry() ;
  UShort_t tIdx = handle.getKey(&tree) ;
  if ( tIdx == FileUtil::NOKEY )
   { tIdx = handle.addTree(tree) ; }
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
// The print format is treeIndex:entryIndex
void CelEventEntry::printEventInfo( const char * /* options */ ) const
 {
  std::cout << _treeIndex << ':' << _entryIndex ;
 }


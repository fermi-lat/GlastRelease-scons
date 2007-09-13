
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

#include "rootUtil/CelEntryIndex.h"
#include "rootUtil/FileUtil.h"
#include "rootUtil/DataHandle.h"
#include "rootUtil/CelFileTreeNames.h"

#include <TTree.h>

#include <iostream>

ClassImp(CelEntryIndex) ;

// Default c'tor. Component Name is not set
CelEntryIndex::CelEntryIndex()
 : BranchGroup(),
   _componentName(),
   _entryIndex(-1,*this,"EntryIndex"),  
   _treeIndex(FileUtil::NOKEY,*this,"TreeIndex")
 { DataHandleInstance::init() ; }
	      
// Standard c'tor.  Component Name is set
CelEntryIndex::CelEntryIndex( const std::string & componentName )
 : BranchGroup(),
   _componentName(componentName),
   _entryIndex(-1,*this,"EntryIndex"),  
   _treeIndex(FileUtil::NOKEY,*this,"TreeIndex")
 { DataHandleInstance::init() ; }	      
	      
// Copy c'tor.  CelComponent Name is copied
CelEntryIndex::CelEntryIndex( const CelEntryIndex & other)
 : BranchGroup(),
   _componentName(other._componentName),
   _entryIndex(other._entryIndex,*this,"EntryIndex"),  
   _treeIndex(other._treeIndex,*this,"TreeIndex")
 { DataHandleInstance::init() ; }	      

// D'tor is a no-op
CelEntryIndex::~CelEntryIndex()
 {}

// Assignment operator
CelEntryIndex & CelEntryIndex::operator=( const CelEntryIndex & other )
 {
  if ( this == &other ) return *this ;
  _componentName = other.componentName() ;
  _entryIndex  = other.entryIndex() ;
  _treeIndex = other.treeIndex() ;
  return *this ;
 }

// Grab the current status from 'tree'
//
// The event index is set from tree.GetReadEntry();
// The tree index is looked up in handle
void CelEntryIndex::set( TTree & tree, CelFileTreeNames & handle )
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
Int_t CelEntryIndex::read( const CelFileTreeNames & handle )
 {
  TTree * t = handle.getTree(_treeIndex) ;
  if ( 0 == t ) return -1 ;
  return t->LoadTree(_entryIndex) ;
 }

// Use 'handle' to look up a tree 
//
// The tree lookup uses _treeIndex to get the tree from 'handle'
TTree * CelEntryIndex::getTree( const CelFileTreeNames & handle ) const
 {
  return handle.getTree(_treeIndex) ;
 }

// Print the information about the current event to cout
// 
// The print format is treeIndex:entryIndex
void CelEntryIndex::printEventInfo( const char * /* options */ ) const
 {
  std::cout << _treeIndex << ':' << _entryIndex ;
 }


// -*- Mode: c++ -*-
#ifndef EventPointer_cxx
#define EventPointer_cxx
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
*
*
*/

//
// EventPointer
//
// EventPointer stores information needed to point to a part of an event
// that is located in another TTree.
//
// EventPointer uses TreeRefHandle to associated TTrees with UShort_t keys.
//
// EventPointer stores two pieces of information for each event.
//
//  _eventIndex (goes to branch XXX_EvtIndex)  is the index of the event component in the TTree
//  _treeKey    (goes to branch XXX_TreeIndex) is the index of the TTree in the associated TreeRefHandle
//

// This Class's header
#include "rootUtil/EventPointer.h"

// c++/stl headers
#include <iostream>

// ROOT Headers
#include <TTree.h>

// Other headers from this package
#include "rootUtil/FileUtil.h"
#include "rootUtil/DataHandle.h"
#include "rootUtil/TreeRefHandle.h"


ClassImp(EventPointer);

EventPointer::EventPointer():
  BranchGroup(),
  _componentName(),
  _eventIndex(-1,*this,"EvtIndex"),  
  _treeKey(FileUtil::NOKEY,*this,"TreeIndex"){
  // Default c'tor.  CelComponent Name is not set
  DataHandleInstance::init();
}
	      
EventPointer::EventPointer(const std::string& componentName):
  BranchGroup(),
  _componentName(componentName),
  _eventIndex(-1,*this,"EvtIndex"),  
  _treeKey(FileUtil::NOKEY,*this,"TreeIndex"){  
  // Standard c'tor.  CelComponent Name is set
  DataHandleInstance::init();
}	      
	      
EventPointer::EventPointer(const EventPointer& other):
  BranchGroup(),
  _componentName(other._componentName),
  _eventIndex(other._eventIndex,*this,"EvtIndex"),  
  _treeKey(other._treeKey,*this,"TreeIndex"){	      
  // Copy c'tor.  CelComponent Name is copied
  DataHandleInstance::init();
}	      

EventPointer::~EventPointer(){
  // D'tor is a no-op
  ;
}

EventPointer& EventPointer::operator=(const EventPointer& other) {
  // Assignment operator
  if ( this == &other ) return *this;
  _componentName = other.componentName();
  _eventIndex  = other.eventIndex();
  _treeKey = other.treeKey();
  return *this;
}

void EventPointer::set(TTree& tree, TreeRefHandle& handle) {
  // Grab the current status from 'tree'
  //
  // The event index is set from tree.GetReadEntry();
  // The tree index is looked up in handle
  _eventIndex = tree.GetReadEntry();
  UShort_t tIdx = handle.getKey(&tree);
  if ( tIdx == FileUtil::NOKEY ) {
    tIdx = handle.addTree(tree);
  }
  _treeKey = tIdx;
}

Int_t EventPointer::read(const TreeRefHandle& handle) {
  // Read an event using handle to look up the correct tree
  //
  // The tree lookup uses _treeKey to get the tree from 'handle'
  // Then the Entry _eventIndex is loaded using LoadTree(_eventIndex);
  TTree* t = handle.getTree(_treeKey);
  if ( 0 == t ) return -1;
  return t->LoadTree(_eventIndex);
}

TTree* EventPointer::getTree(const TreeRefHandle& handle) const {
  // Use 'handle' to look up a tree 
  //
  // The tree lookup uses _treeKey to get the tree from 'handle'
  return handle.getTree(_treeKey);
}

void EventPointer::printEventInfo(const char* /* options */) const {
  // Print the information about this event to cout
  // 
  // The print format is treeKey:eventIndex
  std::cout << _treeKey << ':' <<_eventIndex;
}

#endif

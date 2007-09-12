// -*- Mode: c++ -*-
#ifndef PointerIndex_cxx
#define PointerIndex_cxx
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
// PointerIndex
//
//
// PointerIndex is a sub-class of TVirtualIndex.  It is uses a PointerSkim as an index in to the component trees.
//
// The whole point of this class is to override the method
//   Int_t TVirtualIndex::GetEntryNumberFriend(const TTree* tree);  
//
// Which is called in TTree::LoadTreeFriend(Long64_t entry, TTree* masterTree) to load the friend tree.
// 
// PointerIndex::GetEntryNumberFriend(const TTree* tree) goes into the master tree and finds the correct event index
// in the TChain of the component tree.  
//
// The Index in the component tree is
//
//     GlobalOffset + MetaEventOffset + LocalOffset + EventIndex
//
//     GlobalOffset is the offset from starting in the middle of the chain.  Except in special cases it is 0.
//     MetaEventOffset is the offset from all the previous entries in the FileTree
//     LocalOffset is the offset form all the previous TTrees in the current entry in the FileTree
//     EventIndex is the Index of the event in the tree that it actually lives on.
// 

// This Class's header
#include "rootUtil/PointerIndex.h"

// c++/stl headers
#include <iostream>

// ROOT Headers
#include <TTree.h>

// Other headers from this package
#include "rootUtil/PointerSkim.h"
#include "rootUtil/Component.h"


ClassImp(PointerIndex);

PointerIndex* PointerIndex::buildIndex(PointerSkim& skim, const std::string& component,TTree* tree,Long64_t offset){
  // Build and return pointer index from a pointer skim
  PointerIndex* pIdx = new PointerIndex(skim,component,tree,offset);
  tree->SetTreeIndex(pIdx);
  return pIdx;
}

PointerIndex::PointerIndex():
  TVirtualIndex(),
  _offset(0),
  _skim(0),
  _component(0){
  // Default c'tor.  Does nothing
  SetTree(0);
}

PointerIndex::PointerIndex(PointerSkim& skim, const std::string& compontent, TTree* tree, Long64_t offset):
  TVirtualIndex(),
  _offset(offset),
  _skim(&skim),
  _component(skim.getComponent(skim.componentIndex(compontent))){
  // Standard c'tor.
  SetTree(tree);
}

PointerIndex::~PointerIndex(){
  // Delete stuff
}


Int_t PointerIndex::GetEntryNumberFriend(const TTree* tree){
  // return the index into the component tree for the current event in the master tree

  // Sanity check
  if ( 0 == _skim || 0 == _component ) return -1;

  // Read master tree 
  Long64_t check = _skim->readEventAddress(tree->GetReadEntry());

  // Get the Event index
  Long64_t evtIdx = _offset + _skim->metaOffset();
  evtIdx += _component->getIndexInLocalChain();
  // std::cout << "PointerIndex::GetEntryNumberFriend(" << _skim->eventIndex() << ':' << evtIdx << ')' << std::endl;
  return evtIdx;  
}

Long64_t PointerIndex::GetN() const {
  // Get the total number of Events in the index
  if ( 0 == _skim || 0 == _component ) return 0;
  return _skim->entries();
}

#endif

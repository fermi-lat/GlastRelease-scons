// -*- Mode: c++ -*-
#ifndef CelIndex_cxx
#define CelIndex_cxx
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
// CelIndex
//
//
// CelIndex is a sub-class of TVirtualIndex.  It is uses a CompositeEventList as an index in to the component trees.
//
// The whole point of this class is to override the method
//   Int_t TVirtualIndex::GetEntryNumberFriend(const TTree* tree);  
//
// Which is called in TTree::LoadTreeFriend(Long64_t entry, TTree* masterTree) to load the friend tree.
// 
// CelIndex::GetEntryNumberFriend(const TTree* tree) goes into the master tree and finds the correct event index
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
#include "rootUtil/CelIndex.h"

// c++/stl headers
#include <iostream>

// ROOT Headers
#include <TTree.h>

// Other headers from this package
#include "rootUtil/CompositeEventList.h"
#include "rootUtil/CelComponent.h"


ClassImp(CelIndex);

CelIndex* CelIndex::buildIndex(CompositeEventList& skim, const std::string& component,TTree* tree,Long64_t offset){
  // Build and return an index from a composite event list
  CelIndex* pIdx = new CelIndex(skim,component,tree,offset);
  tree->SetTreeIndex(pIdx);
  return pIdx;
}

CelIndex::CelIndex():
  TVirtualIndex(),
  _offset(0),
  _skim(0),
  _component(0){
  // Default c'tor.  Does nothing
  SetTree(0);
}

CelIndex::CelIndex(CompositeEventList& skim, const std::string& compontent, TTree* tree, Long64_t offset):
  TVirtualIndex(),
  _offset(offset),
  _skim(&skim),
  _component(skim.getComponent(skim.componentIndex(compontent))){
  // Standard c'tor.
  SetTree(tree);
}

CelIndex::~CelIndex(){
  // Delete stuff
}


Int_t CelIndex::GetEntryNumberFriend(const TTree* tree){
  // return the index into the component tree for the current event in the master tree

  // Sanity check
  if ( 0 == _skim || 0 == _component ) return -1;

  // Read master tree 
  Long64_t check ;
  check = _skim->readEventAddress(tree->GetReadEntry());

  // Get the Event index
  Long64_t evtIdx = _offset + _skim->setOffset();
  evtIdx += _component->getIndexInLocalChain();
  // std::cout << "CelIndex::GetEntryNumberFriend(" << _skim->eventIndex() << ':' << evtIdx << ')' << std::endl;
  return evtIdx;  
}

Long64_t CelIndex::GetN() const {
  // Get the total number of Events in the index
  if ( 0 == _skim || 0 == _component ) return 0;
  return _skim->entries();
}

#endif

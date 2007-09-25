
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

#include "rootUtil/CelIndex.h"
#include "rootUtil/CompositeEventList.h"
#include "rootUtil/CelEventComponent.h"
#include <TTree.h>
#include <Riostream.h>

ClassImp(CelIndex) ;

CelIndex * CelIndex::buildIndex
 ( CompositeEventList & cel,
   const TString & componentName, TTree * tree, Long64_t offset )
 {
  CelIndex * pIdx = new CelIndex(cel,componentName,tree,offset) ;
  tree->SetTreeIndex(pIdx) ;
  return pIdx ;
 }

CelIndex::CelIndex()
 : TVirtualIndex(),
   _offset(0), _cel(0), _component(0)
 { SetTree(0) ; }

CelIndex::CelIndex
 ( CompositeEventList & cel,
   const TString & componentName,
   TTree * tree, Long64_t offset )
 : TVirtualIndex(),
   _offset(offset), _cel(&cel),
   _component(cel.getComponent(cel.componentIndex(componentName)))
 { SetTree(tree) ; }

CelIndex::~CelIndex()
 {}


Int_t CelIndex::GetEntryNumberFriend(const TTree* tree){
  // return the index into the component tree for the current event in the master tree

  // Sanity check
  if ( 0 == _cel || 0 == _component ) return -1;

  // Read master tree 
  Long64_t check ;
  check = _cel->readEventAddress(tree->GetReadEntry());

  // Get the Event index
  Long64_t evtIdx = _offset + _cel->setOffset();
  evtIdx += _component->getIndexInLocalChain();
  // std::cout << "CelIndex::GetEntryNumberFriend(" << _cel->eventIndex() << ':' << evtIdx << ')' << std::endl;
  return evtIdx;  
}

Long64_t CelIndex::GetN() const {
  // Get the total number of Events in the index
  if ( 0 == _cel || 0 == _component ) return 0;
  return _cel->nbEvents();
}



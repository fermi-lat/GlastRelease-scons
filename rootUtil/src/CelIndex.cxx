
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
   _component(cel.getComponent(componentName))
 { SetTree(tree) ; }

CelIndex::~CelIndex()
 {}


Int_t CelIndex::GetEntryNumberFriend( const TTree * tree )
 {
  // return the index into the component tree for the current event in the master tree

  // Sanity check
  if (_cel==0||_component==0) { return -1 ; }

  // Read entry in cel trees
  Long64_t check ;
  check = _cel->shallowRead(tree->GetReadEntry()) ;

  // Get the Event index
  Long64_t evtIdx = _offset + _cel->fileSetOffset() ;
  evtIdx += _component->getIndexInLocalChain() ;
  // std::cout << "CelIndex::GetEntryNumberFriend(" << _cel->currentEventIndex() << ':' << evtIdx << ')' << std::endl;
  return evtIdx ;  
 }

Long64_t CelIndex::GetN() const
 {
  if (_cel==0||_component==0) { return 0 ; }
  return _cel->numEvents() ;
 }



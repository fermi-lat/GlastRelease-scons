// -*- Mode: c++ -*-
#ifndef CelEventComponent_cxx
#define CelEventComponent_cxx
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
// CelEventComponent stores all the information needed to point to a part of an event
// that is located in another TTree.
//
// This for a single event this information is stored on two TTees. 
//   The event information: index and tree key is stored on one tree
//   The collection information: lists of trees is stored on the other tree
//
// The information to navigate between these two tree is stored by the CompositeEventList class
// since in is shared by all the components in a skim.
//


// This Class's header
#include "rootUtil/CelEventComponent.h"

// c++/stl headers
#include <iostream>

// ROOT Headers
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TEventList.h>

// Other headers from this package
#include "rootUtil/OptUtil.h"
#include "rootUtil/DataHandle.h"
#include "rootUtil/CelFileAndTreeSet.h"

ClassImp(CelEventComponent);


CelEventComponent::CelEventComponent():
  _componentName("NULL"),
  _event(),
  _tree(){
  // Default c'tor.  Needed for ROOT
}
	      
CelEventComponent::CelEventComponent(const std::string& componentName):
  _componentName(componentName),
  _event(componentName),
  _tree(componentName){
  // Standard c'tor, stores the name of the component 
}	      

CelEventComponent::~CelEventComponent(){
  // D'tor
  ;
}


void CelEventComponent::set(TTree& tree) {
  // set the current event
  _event.set(tree,_tree);
}

Int_t CelEventComponent::read() {
  // read an event
  return _event.read(_tree);
}

TTree* CelEventComponent::getTree() const {
  // Get the Tree that is being read
  return _event.getTree(_tree);
}


Int_t CelEventComponent::makeBranches( TTree & fileTree, TTree & eventTree, Int_t bufsize) const {
  Int_t n_e = _event.makeBranches(eventTree,_componentName.c_str(),bufsize);
  if ( n_e < 0 ) return n_e;
  Int_t n_f = _tree.makeBranches(fileTree,_componentName.c_str(),bufsize);
  if ( n_f < 0 ) return n_f;
  return n_e + n_f;
}

Int_t CelEventComponent::attachToTree(TTree& fileTree, TTree& eventTree ) {
  Int_t n_e = _event.attachToTree(eventTree,_componentName.c_str());
  if ( n_e < 0 ) return n_e;
  Int_t n_f = _tree.attachToTree(fileTree,_componentName.c_str());
  if ( n_f < 0 ) return n_f;
  return n_e + n_f; 
}


Bool_t CelEventComponent::addToChain(TChain*& chain) {
  // Building a TChain
  for ( UShort_t iT(0); iT < _tree.size(); iT++ ) {
    TTree* t = _tree.getTree(iT);
    if ( 0 == t ) {
      return kFALSE;
    }
    if ( chain == 0 ) {
      chain = new TChain(t->GetName(),"Working on the Chain Gain");      
    }
    chain->Add(t->GetCurrentFile()->GetName());
  }
  return kTRUE;
} 


Long64_t CelEventComponent::getLocalOffset() const {
  // Get the offset to the first event
  Long64_t localOffset = _tree.getOffset(_event.treeIndex());
  return localOffset;
}

Long64_t CelEventComponent::getIndexInLocalChain() const {
  // Get the index of the event the local chain (ie, in this meta event)
  Long64_t evtIdx = getLocalOffset();
  evtIdx += _event.entryIndex();
  return evtIdx;
}


void CelEventComponent::dumpEvent(const char* options) const {
  if ( OptUtil::has_option(options,'v') ) {
    _tree.printTreeInfo(_event.treeIndex(),options);
  }  
  _event.printEventInfo(options);
}

void CelEventComponent::listTrees(const char* options) const {
  _tree.show(options);
}

#endif

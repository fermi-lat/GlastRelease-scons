
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

#include <rootUtil/CelEventComponent.h>
#include <rootUtil/OptUtil.h>
#include <rootUtil/DataHandle.h>
#include <rootUtil/CelFileAndTreeSet.h>

#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TEventList.h>
#include <Riostream.h>



//====================================================================
// 
//====================================================================


ClassImp(CelEventComponent) ;

CelEventComponent::CelEventComponent()
 : _componentName("NULL"), _currentEntryIndex(), _currentSet()
 {}
	      
CelEventComponent::CelEventComponent( const TString & componentName )
 : _componentName(componentName), _currentEntryIndex(componentName), _currentSet(componentName)
 {}	      

CelEventComponent::~CelEventComponent()
 {}



//====================================================================
// 
//====================================================================


void CelEventComponent::set( TTree& tree) {
  // set the current event
  _currentEntryIndex.set(tree,_currentSet);
}

Int_t CelEventComponent::read() {
  // read an event
  return _currentEntryIndex.read(_currentSet);
}

TTree* CelEventComponent::getTree() const {
  // Get the Tree that is being read
  return _currentEntryIndex.getTree(_currentSet);
}


Int_t CelEventComponent::makeBranches( TTree & fileTree, TTree & eventTree, Int_t bufsize) const {
  Int_t n_e = _currentEntryIndex.makeBranches(eventTree,_componentName.Data(),bufsize);
  if ( n_e < 0 ) return n_e;
  Int_t n_f = _currentSet.makeBranches(fileTree,_componentName.Data(),bufsize);
  if ( n_f < 0 ) return n_f;
  return n_e + n_f;
}

Int_t CelEventComponent::attachToTree(TTree& fileTree, TTree& eventTree ) {
  Int_t n_e = _currentEntryIndex.attachToTree(eventTree,_componentName.Data());
  if ( n_e < 0 ) return n_e;
  Int_t n_f = _currentSet.attachToTree(fileTree,_componentName.Data());
  if ( n_f < 0 ) return n_f;
  return n_e + n_f; 
}

// Building a TChain
Bool_t CelEventComponent::addToChain( TChain * & chain )
 {
  // TO REVIEW WITH ERIC : DO WE REALLY NEED TO GET THE REAL TREE
  // RATHER THAN PLAY ONLY WITH NAMES ??
  UShort_t iT ;
  for ( iT=0 ; iT<_currentSet.size() ; iT++ )
   {
    TTree * t = _currentSet.getTree(iT) ;
    if (t==0) { return kFALSE ; }
    if (chain==0)
     { chain = new TChain(t->GetName(),"Cel Component Chain") ; }
    chain->Add(t->GetCurrentFile()->GetName()) ;
   }
  return kTRUE ;
 } 


Long64_t CelEventComponent::getLocalOffset() const {
  // Get the offset to the first event
  Long64_t localOffset = _currentSet.getOffset(_currentEntryIndex.treeIndex());
  return localOffset;
}

Long64_t CelEventComponent::getIndexInLocalChain() const {
  // Get the index of the event the local chain (ie, in this meta event)
  Long64_t evtIdx = getLocalOffset();
  evtIdx += _currentEntryIndex.entryIndex();
  return evtIdx;
}

void CelEventComponent::dumpEvent( const char * options ) const
 {
  if ( OptUtil::has_option(options,'v') )
   { _currentSet.printTreeInfo(_currentEntryIndex.treeIndex(),options) ; }  
  _currentEntryIndex.printEventInfo(options) ;
 }

void CelEventComponent::listTrees( const char * options ) const
 { _currentSet.show(options) ; }



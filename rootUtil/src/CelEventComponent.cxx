
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
#include <rootUtil/BgDataHandle.h>
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


void CelEventComponent::set( TTree & tree)
 {
  // set the current event
  _currentEntryIndex.set(tree,_currentSet) ;
 }

Int_t CelEventComponent::read()
 {
  // read an event
  return _currentEntryIndex.read(_currentSet) ;
 }

TTree * CelEventComponent::getTree() const
 {
  // Get the Tree that is being read
  return _currentEntryIndex.getTree(_currentSet) ;
 }


Int_t CelEventComponent::makeBranches( TTree & fileTree, TTree & eventTree, Int_t bufsize) const
 {
  Int_t n_e = _currentEntryIndex.makeBranches(eventTree,_componentName.Data(),bufsize);
  if ( n_e < 0 ) return n_e ;
  Int_t n_f = _currentSet.makeBranches(fileTree,_componentName.Data(),bufsize);
  if ( n_f < 0 ) return n_f ;
  return n_e + n_f ;
 }

Int_t CelEventComponent::attachToTree( TTree & entryTree, TTree & fileTree )
 {
  Int_t n_e = _currentEntryIndex.attachToTree(entryTree,_componentName.Data()) ;
  if ( n_e < 0 ) return n_e ;
  Int_t n_f = _currentSet.attachToTree(fileTree,_componentName.Data()) ;
  if ( n_f < 0 ) return n_f ;
  return n_e + n_f; 
 }

// Building a TChain
Bool_t CelEventComponent::addToChain( TChain * & chain )
 { return _currentSet.addToChain(chain) ; }

Long64_t CelEventComponent::indexInCurrentSet() const
 {
  UShort_t treeIndex = _currentEntryIndex.treeIndex() ;
  Long64_t previousTreesOffset = _currentSet.getOffset(treeIndex) ;
  Long64_t evtIndex = previousTreesOffset + _currentEntryIndex.entryIndex() ;
  return evtIndex ;
 }

void CelEventComponent::dumpEvent( const char * options ) const
 {
  if ( OptUtil::has_option(options,'v') )
   { _currentSet.printTreeInfo(_currentEntryIndex.treeIndex(),options) ; }  
  _currentEntryIndex.printEventInfo(options) ;
 }

void CelEventComponent::listTrees( const char * options ) const
 { _currentSet.show(options) ; }



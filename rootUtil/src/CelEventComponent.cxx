
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
 : _componentName("NULL"), _currentEntryIndex(), _currentSet(),
   _tree(0), _data(0)
 {}
	      
CelEventComponent::CelEventComponent( const TString & componentName )
 : _componentName(componentName), _currentEntryIndex(componentName), _currentSet(componentName)
 {}	      

CelEventComponent::~CelEventComponent()
 {}



//====================================================================
// WRITING
//====================================================================


void CelEventComponent::registerEntry( TTree & tree )
 { _currentEntryIndex.set(tree,_currentSet) ; }

void CelEventComponent::nextSet()
 {
  _currentOffset.increment(_currentSet.entries()) ;
  _currentSet.reset() ;
 }



//====================================================================
// READING
//====================================================================


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


Int_t CelEventComponent::makeBranches( TTree * entryTree, TTree * fileTree, TTree * offsetTree, Int_t bufsize) const
 {
  Int_t n, nTot = 0 ;
  if (entryTree!=0)
   {
    n = _currentEntryIndex.makeBranches(*entryTree,_componentName.Data(),bufsize) ;
    if (n<0) return n ;
    nTot += n ;
   }
  if (fileTree!=0)
   {
	n = _currentSet.makeBranches(*fileTree,_componentName.Data(),bufsize) ;
	if (n<0) return n ;
	nTot += n ;
   }
  if (offsetTree!=0)
   {
	n = _currentOffset.makeBranches(*offsetTree,_componentName.Data(),bufsize) ;
	if (n<0) return n ;
	nTot += n ;
   }
  return nTot ; 
 }

Int_t CelEventComponent::attachToTree( TTree * entryTree, TTree * fileTree, TTree * offsetTree )
 {
  Int_t n, nTot = 0 ;
  if (entryTree!=0)
   {
	n = _currentEntryIndex.attachToTree(*entryTree,_componentName.Data()) ;
	if (n<0) return n ;
	nTot += n ;
   }
  if (fileTree!=0)
   {
	n = _currentSet.attachToTree(*fileTree,_componentName.Data()) ;
	if (n<0) return n ;
	nTot += n ;
   }
  if (offsetTree!=0)
   {
	n = _currentOffset.attachToTree(*offsetTree,_componentName.Data()) ;
	if (n<0) return n ;
	nTot += n ;
   }
  return nTot ; 
 }

// Building a TChain
Bool_t CelEventComponent::addToChain( TChain * & chain )
 { return _currentSet.addToChain(chain) ; }

Long64_t CelEventComponent::currentIndexInChain() const
 {
  UShort_t treeIndex = _currentEntryIndex.treeIndex() ;
  Long64_t previousTreesOffset = _currentSet.getOffset(treeIndex) ;
  Long64_t previousSetsOffset = _currentOffset.getOffset() ;
  Long64_t evtIndex = previousSetsOffset + previousTreesOffset + _currentEntryIndex.entryIndex() ;
  //std::cout<<"("<<previousSetsOffset<<"+"<<previousTreesOffset<<"+"<<_currentEntryIndex.entryIndex()<<") "<<std::flush ;
  return evtIndex ;
 }

void CelEventComponent::printEventInfo( const char * options ) const
 {
  if ( OptUtil::has_option(options,'v') )
   { _currentSet.printTreeInfo(_currentEntryIndex.treeIndex(),options) ; }  
  _currentEntryIndex.printInfo() ;
 }

void CelEventComponent::printSetInfo( const char * options, const char * prefix ) const
 { _currentSet.printTreesInfo("tf",prefix) ; }



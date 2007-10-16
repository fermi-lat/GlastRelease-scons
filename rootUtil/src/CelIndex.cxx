
/*
* Project: GLAST
* Package: rootUtil
*    File: $Id$
* Authors:
*   EC, Eric Charles, SLAC, echarles@slac.stanford.edu
*   DC, David Chamont, LLR, chamont@llr.in2p3.fr
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

CelIndex::CelIndex()
 : TVirtualIndex(), _cel(0), _component(0)
 { SetTree(0) ; }

CelIndex::CelIndex
 ( CompositeEventList & cel,
   const TString & componentName,
   TChain * componentChain )
 : TVirtualIndex(), _cel(&cel),
   _component(cel.getComponent(componentName))
 { SetTree((TTree*)componentChain) ; }

CelIndex::~CelIndex()
 {}

Long64_t CelIndex::GetN() const
 {
  if (_cel==0||_component==0) { return 0 ; }
  return _cel->numEvents() ;
 }

Int_t CelIndex::GetEntryNumberFriend( const TTree * celCompositeEvents )
 {
  if (_cel==0||_component==0) { return -1 ; }
  _cel->shallowRead(celCompositeEvents->GetReadEntry()) ;
  Long64_t evtIdx = _cel->fileSetOffset() ;
  evtIdx += _component->indexInCurrentSet() ;
  // std::cout << "CelIndex::GetEntryNumberFriend(" << _cel->currentEventIndex() << '|' << evtIdx << ')' << std::endl;
  return evtIdx ;  
 }



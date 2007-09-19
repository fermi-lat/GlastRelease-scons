
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

#include "rootUtil/CelLink.h"
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TObjArray.h>
#include <iostream>
#include <assert.h>

ClassImp(CelLink) ;

CelLink::CelLink()
 : BranchGroup(),
   _eventIndex(-1,*this,"EventIndex"),
   _setIndex(0,*this,"SetIndex"),
   _setOffset(0,*this,"SetOffset")
 { DataHandleInstance::init() ; }

CelLink::CelLink( const CelLink & other )
 : BranchGroup(),
   _eventIndex(other._eventIndex,*this,"EventIndex"),
   _setIndex(other._setIndex,*this,"SetIndex"),
   _setOffset(other._setOffset,*this,"SetOffset")
 { DataHandleInstance::init() ; }	      

CelLink::~CelLink()
 {}

CelLink & CelLink::operator=( const CelLink & other )
 {
  if ( this == &other ) return *this ;
  _eventIndex = other._eventIndex ;
  _setIndex  = other._setIndex ;
  _setOffset = other._setOffset ;
  return *this ;
 }

// Print the information about the current event to cout
// 
// The print format is eventIndex:setIndex:setOffset
void CelLink::printEventInfo( const char * /* options */ ) const
 {
  std::cout << _eventIndex << ':' << _setIndex << ':' << _setOffset ;
 }


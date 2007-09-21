
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

#include "rootUtil/CelEventLink.h"
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TObjArray.h>
#include <iostream>
#include <assert.h>

ClassImp(CelEventLink) ;

CelEventLink::CelEventLink()
 : BranchGroup(),
   _eventIndex(-1,*this,"EventIndex"),
   _setIndex(0,*this,"SetIndex"),
   _setOffset(0,*this,"SetOffset")
 { DataHandleInstance::init() ; }

CelEventLink::CelEventLink( const CelEventLink & other )
 : BranchGroup(),
   _eventIndex(other._eventIndex,*this,"EventIndex"),
   _setIndex(other._setIndex,*this,"SetIndex"),
   _setOffset(other._setOffset,*this,"SetOffset")
 { DataHandleInstance::init() ; }	      

CelEventLink::~CelEventLink()
 {}

CelEventLink & CelEventLink::operator=( const CelEventLink & other )
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
void CelEventLink::printEventInfo( const char * /* options */ ) const
 {
  std::cout << _eventIndex << ':' << _setIndex << ':' << _setOffset ;
 }


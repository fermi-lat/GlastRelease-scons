
/*
* Project: GLAST
* Package: rootUtil
*    File: $Id$
* Authors:
*   EC, Eric Charles,    SLAC              echarles@slac.stanford.edu
*   DC, David Chamont, LLR, chamont@llr.in2p3.fr
*
* Copyright (c) 2007
*   Regents of Stanford University. All rights reserved.
*
*/

#include <rootUtil/CelEventLink.h>
#include <Riostream.h>
#include <assert.h>

ClassImp(CelEventLink) ;

CelEventLink::CelEventLink()
 : BranchGroup(),
   _eventIndex(-1,*this,"Event_Index"),
   _setIndex(0,*this,"Event_SetIndex")
 { BgDataHandleInstance::init() ; }

CelEventLink::CelEventLink( const CelEventLink & other )
 : BranchGroup(),
   _eventIndex(other._eventIndex,*this,"Event_Index"),
   _setIndex(other._setIndex,*this,"Event_SetIndex")
 {}	      

CelEventLink::~CelEventLink()
 {}

CelEventLink & CelEventLink::operator=( const CelEventLink & other )
 {
  if ( this == &other ) return *this ;
  _eventIndex = other._eventIndex ;
  _setIndex  = other._setIndex ;
  return *this ;
 }

void CelEventLink::printInfo( const char * /* options */ ) const
 { std::cout << _eventIndex << '|' << _setIndex ; }


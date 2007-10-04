
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
   _fileSetIndex(0,*this,"Event_SetIndex"),
   _fileSetOffset(0,*this,"Event_SetOffset")
 { BgDataHandleInstance::init() ; }

CelEventLink::CelEventLink( const CelEventLink & other )
 : BranchGroup(),
   _eventIndex(other._eventIndex,*this,"Event_Index"),
   _fileSetIndex(other._fileSetIndex,*this,"Event_SetIndex"),
   _fileSetOffset(other._fileSetOffset,*this,"Event_SetOffset")
 { BgDataHandleInstance::init() ; }	      

CelEventLink::~CelEventLink()
 {}

CelEventLink & CelEventLink::operator=( const CelEventLink & other )
 {
  if ( this == &other ) return *this ;
  _eventIndex = other._eventIndex ;
  _fileSetIndex  = other._fileSetIndex ;
  _fileSetOffset = other._fileSetOffset ;
  return *this ;
 }

void CelEventLink::printInfo( const char * /* options */ ) const
 {
  std::cout
    << _eventIndex
    << ':' << _fileSetIndex
    << ':' << _fileSetOffset ;
 }


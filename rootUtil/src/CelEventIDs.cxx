
/*
* Project: GLAST
* Package: rootUtil
*    File: $Id$
* Authors:
*   DC, David Chamont, LLR, chamont@llr.in2p3.fr
*
* Copyright (c) 2007
*                   Regents of Stanford University. All rights reserved.
*
*
*
*/

#include <rootUtil/CelEventIDs.h>
#include <Riostream.h>
#include <cassert>
#include <cstring>

ClassImp(CelEventIDs) ;

CelEventIDs::CelEventIDs()
 : BranchGroup(),
   _productionID(0,*this,"Event_ProductionID"),
   _runID(-1,*this,"Event_RunID"),
   _eventID(-1,*this,"Event_EventID")
 { BgDataHandleInstance::init() ; }

CelEventIDs::CelEventIDs( const CelEventIDs & other )
 : BranchGroup(),
   _productionID(strdup(other._productionID),*this,"Event_ProductionID"),
   _runID(other._runID,*this,"Event_RunID"),
   _eventID(other._eventID,*this,"Event_EventID")
 { BgDataHandleInstance::init() ; }	      

CelEventIDs::~CelEventIDs()
 {}

CelEventIDs & CelEventIDs::operator=( const CelEventIDs & other )
 {
  if ( this == &other ) return *this ;
  delete _productionID ;
  _productionID = strdup(other._productionID) ;
  _runID = other._runID ;
  _eventID = other._eventID ;
  return *this ;
 }

void CelEventIDs::printInfo( const char * /* options */ ) const
 {
  std::cout
    << _productionID
    << '|' << _runID
    << '|' << _eventID ;
 }


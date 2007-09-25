
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

#include "rootUtil/CelEventInfo.h"
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TObjArray.h>
#include <iostream>
#include <assert.h>

ClassImp(CelEventInfo) ;

CelEventInfo::CelEventInfo()
 : BranchGroup(),
   _eventIndex(-1,*this,"Index"),
   _runId(-1,*this,"RunId"),
   _eventId(-1,*this,"EventId"),
   _fileSetIndex(0,*this,"FileSetIndex"),
   _fileSetOffset(0,*this,"FileSetOffset")
 { DataHandleInstance::init() ; }

CelEventInfo::CelEventInfo( const CelEventInfo & other )
 : BranchGroup(),
   _eventIndex(other._eventIndex,*this,"Index"),
   _runId(other._runId,*this,"RunId"),
   _eventId(other._eventId,*this,"EventId"),
   _fileSetIndex(other._fileSetIndex,*this,"FileSetIndex"),
   _fileSetOffset(other._fileSetOffset,*this,"FileSetOffset")
 { DataHandleInstance::init() ; }	      

CelEventInfo::~CelEventInfo()
 {}

CelEventInfo & CelEventInfo::operator=( const CelEventInfo & other )
 {
  if ( this == &other ) return *this ;
  _eventIndex = other._eventIndex ;
  _runId = other._runId ;
  _eventId = other._eventId ;
  _fileSetIndex  = other._fileSetIndex ;
  _fileSetOffset = other._fileSetOffset ;
  return *this ;
 }

void CelEventInfo::printEventInfo( const char * /* options */ ) const
 {
  std::cout
    << _eventIndex
    << ':' << _runId
    << ':' << _eventId
    << ':' << _fileSetIndex
    << ':' << _fileSetOffset ;
 }


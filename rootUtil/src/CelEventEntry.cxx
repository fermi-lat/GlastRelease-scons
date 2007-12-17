
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

#include <rootUtil/CelEventEntry.h>
#include <rootUtil/FileUtil.h>
#include <rootUtil/BgDataHandle.h>
#include <rootUtil/CelFileAndTreeSet.h>

#include <TTree.h>
#include <Riostream.h>



//====================================================================
// Construction and assignment
// One of the use-case is that one could want to manipulate
// two CELs at the same time, one for reading and one for
// writing. In this context, it could make sense to construct
// first the reading CEL, and duplicate it so to get a
// write one with same structure. This is why we have copy
// constructors which are copying the BgDataHandle values
// but relinking them to the new BranchGroup. Also, the assignment
// is assigning only the value, not the link between handles and group,
// so that ultimately we could do writeCel = readCel which would
// duplicate values of the current Event.
// As a summary : copy construction is duplicating the structure,
// and assigment only concerns the current event values.
//====================================================================


ClassImp(CelEventEntry) ;

CelEventEntry::CelEventEntry()
 : _componentName(),
   _treeIndex(FileUtil::NOKEY,*this,"Event_TreeIndex"),
   _entryIndex(-1,*this,"Event_EntryIndex")
 { BgDataHandleInstance::init() ; }
	      
CelEventEntry::CelEventEntry( const TString & componentName )
 : _componentName(componentName),
   _treeIndex(FileUtil::NOKEY,*this,"Event_TreeIndex"),
   _entryIndex(-1,*this,"Event_EntryIndex")
 { BgDataHandleInstance::init() ; }	      
	      
CelEventEntry::CelEventEntry( const CelEventEntry & other )
 : BranchGroup(),
   _componentName(other._componentName),
   _treeIndex(other._treeIndex,*this,"Event_TreeIndex"),
   _entryIndex(other._entryIndex,*this,"Event_EntryIndex")
 {}	      

CelEventEntry::~CelEventEntry()
 {}

CelEventEntry & CelEventEntry::operator=( const CelEventEntry & other )
 {
  if ( this == &other ) return *this ;
  _componentName = other.componentName() ; // why not assert it is the same ?
  _treeIndex = other.treeIndex() ;
  _entryIndex  = other.entryIndex() ;
  return *this ;
 }



//====================================================================
// Other
//====================================================================


void CelEventEntry::set( UShort_t treeIndex, Long64_t entryIndex )
 { _treeIndex = treeIndex ;  _entryIndex = entryIndex ; }

void CelEventEntry::printInfo() const
 {
  std::cout
    <<setw(2)<<std::right<< _treeIndex
    <<'|'
    <<setw(2)<<std::left<< _entryIndex ;
 }


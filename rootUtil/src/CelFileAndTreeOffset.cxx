
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


#include "rootUtil/CelFileAndTreeOffset.h"

// c++/stl headers
#include <iostream>
#include <assert.h>

// ROOT Headers
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TObjString.h>
#include <TObjArray.h>
//#include <TArrayL64.h>

// Other headers from this package
#include "rootUtil/RuUtil.h"
#include "rootUtil/BgDataHandle.h"
#include "rootUtil/RuUtil.h"


ClassImp(CelFileAndTreeOffset) ;

CelFileAndTreeOffset::CelFileAndTreeOffset()
 : _componentName(),
   _offset(0,*this,"Set_EntryOffset")
 { BgDataHandleInstance::init() ; }

CelFileAndTreeOffset::CelFileAndTreeOffset( const TString & componentName )
 : _componentName(componentName),
 _offset(0,*this,"Set_EntryOffset")
 { BgDataHandleInstance::init() ; }

CelFileAndTreeOffset::CelFileAndTreeOffset( const CelFileAndTreeOffset & other )
 : BranchGroup(),
   _componentName(other._componentName),
   _offset(other._offset,*this,"Set_EntryOffset")
 {}	      

CelFileAndTreeOffset & CelFileAndTreeOffset::operator=( const CelFileAndTreeOffset & other )
 {
  if ( this == &other ) return *this ;
  _componentName = other.componentName() ; // why not assert it is the same ?
  _offset  = other.getOffset() ;
  return *this ;
 }

CelFileAndTreeOffset::~CelFileAndTreeOffset()
 {}



//======================================================================================
//
// Printing & debugging
//
//======================================================================================


// Print the information about the current set to cout
// 
// The print format is treeIndex|entryIndex
void CelFileAndTreeOffset::printInfo( const char * /* options */ ) const
 {
  std::cout << _offset ;
 }

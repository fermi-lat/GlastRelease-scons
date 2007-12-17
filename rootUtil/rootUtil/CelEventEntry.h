
#ifndef CelEventEntry_h
#define CelEventEntry_h

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

#include "BranchGroup.h"
#include "BgDataHandle.h"
#include <TString.h>
class CelFileAndTreeSet ;

//
// CelEventEntry stores information needed to point to an event component,
// which is an entry located in another tree. It includes the entry index,
// and the index of the tree in the current set of files and trees.
// 

class CelEventEntry : public BranchGroup
 {
	
  public:

    // c'tors & d'tor
    CelEventEntry() ; // Needed for ROOT
    CelEventEntry( const TString & componentName ) ;
    CelEventEntry( const CelEventEntry & ) ;
    CelEventEntry & operator=( const CelEventEntry & ) ;
    ~CelEventEntry() ;

    // Accessors & modifiers
    inline const TString & componentName() const { return _componentName ; }
    inline UShort_t treeIndex() const { return _treeIndex ; }
    inline Long64_t entryIndex() const { return _entryIndex ; }
    void set( UShort_t treeIndex, Long64_t entryIndex ) ;

    // Print the info about the entry
    void printInfo() const ;  

  private :
  
    // Data
	TString _componentName ;
    BgDataHandle<UShort_t> _treeIndex ;  //! Index of the current data TTree in the associated CelFileAndTreeSet
    BgDataHandle<Long64_t> _entryIndex ; //! Index of the current entry in the current data TTree
  
    ClassDef(CelEventEntry,0)

 } ;


#endif

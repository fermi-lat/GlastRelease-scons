
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

    // Access
    inline const TString & componentName() const { return _componentName ; }
    inline Long64_t entryIndex() const { return _entryIndex ; }
    inline UShort_t treeIndex() const { return _treeIndex ; }

    // Methods and functions
    // set the current event
    void set( TTree & tree, CelFileAndTreeSet & handle) ;
    // read an event
    Int_t read( const CelFileAndTreeSet & handle ) ;
    // get the tree that is being read
    TTree * getTree( const CelFileAndTreeSet & handle ) const ;

    // print the info about the entry
    void printInfo() const ;  

  private :
  
    // Data
	TString _componentName ;
    BgDataHandle<Long64_t> _entryIndex ; //! Index of the current entry in the current data TTree
    BgDataHandle<UShort_t> _treeIndex ;  //! Index of the current data TTree in the associated CelFileAndTreeSet
  
    ClassDef(CelEventEntry,0)

 } ;


#endif

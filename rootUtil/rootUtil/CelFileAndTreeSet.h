
#ifndef CelFileAndTreeSet_h
#define CelFileAndTreeSet_h

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
class RuChain ;
class TTree ;
class TFile ;
class TObjArray ;
class TArrayL64 ;
class TObjString ;

#include <map>

//
// CelFileAndTreeSet stores information needed to point to a part of an event
// that is located in another TTree.
//
// CelFileAndTreeSet keep track of the TTrees where the event components live.
//
// Lists of files and tree names are stored persistently on a TTree by 
// CelFileAndTreeSet.  These are:
//
//  _fileNames (goes to branch XXX_FileNames) is the vector of file names
//  _treeNames (goes to branch XXX_TreeNames) is the vector of tree names
// 
// At any given time a CelFileAndTreeSet may be keeping track of any number of TTree
//   
// Access to the TTree is by key
//  TTree* getTree(UShort_t key) const ;
//
// The keys are assigned when TTree are first associated with the CelFileAndTreeSet
// 

class CelFileAndTreeSet  : public BranchGroup
 {

  public :

    // construction
    CelFileAndTreeSet() ; // Needed for ROOT
    CelFileAndTreeSet( const TString & componentName ) ;
    virtual ~CelFileAndTreeSet() ;

    // Methods and functions
    // Reset this object and clear caches
    void resetAll() ;
    void resetCache() ;
    
    // Data Access
    UShort_t addTree( TTree & tree) ;
    UShort_t addTree( const TString & filePath, const TString & treeName, Long64_t nbTreeEntries ) ;
    UShort_t getKey( TTree * tree) const ;
    UShort_t getKey( const TString & filePath, const TString & treeName ) const ;
    TTree * getTree( UShort_t key ) const ;
    const TObjString * getFileName( UShort_t key ) const ;
    const TObjString * getTreeName( UShort_t key ) const ;
    Long64_t getTreeNbEntries( UShort_t key ) const ;
    Long64_t getOffset( UShort_t key ) const ;

    // reading interface
    Bool_t addToChain( RuChain * & chain ) ;  
  
    // Override the methods in BranchGroup.  
    virtual Int_t makeBranches( TTree & celTree, const char* prefix = 0, Int_t bufsize = 32000) const ;
    virtual Int_t attachToTree( TTree & celTree, const char* prefix = 0) ;

    // Access
    inline UShort_t size() const { return _setSize ; }
    inline Long64_t entries() const { return _treesSize ; }
    inline const TString & componentName() const { return _componentName ; }

    // Printing
    // Print the list of trees
    void printTreesInfo( const char* options ="tf", const char* prefix ="" ) const ;
    // Print information about a single tree
    void printTreeInfo( UShort_t key, const char* options ="tf", const char* prefix =""  ) const ;

  private :

    // disable copying and assignment
    CelFileAndTreeSet( const CelFileAndTreeSet & ) ;
    CelFileAndTreeSet & operator=( const CelFileAndTreeSet & ) ;

    // Utility function to actually go and get a tree out of a file
    TTree * fetchTree( UShort_t key ) const ;

    // data
    TString _componentName ;

    // data on branches
    BgDataHandle<UShort_t>  _setSize ;     //! Number of trees stored by this CelFileAndTreeSet
    TObjArray*            _fileNames ;   //! Names of files where the TTree live
    TObjArray*            _treeNames ;   //! Names of trees stored by this CelFileAndTreeSet
    BgDataHandle<Long64_t>  _treesSize ;   //! Total number of entries on the trees stored by this CelFileAndTreeSet
    TArrayL64*            _treeOffsets ; //! Number of entries in TTrees

    // Cache for fast lookup
    mutable std::map<UShort_t,TTree*> _cache ;  //! Maps key -> tree
    mutable std::map<TString,UShort_t> _lookup ; //! Reverse map, for making sure we don't screw up

    ClassDef(CelFileAndTreeSet,0) // For handling referenece to trees
 } ;


#endif


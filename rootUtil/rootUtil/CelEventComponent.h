
#ifndef CelEventComponent_h
#define CelEventComponent_h

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

#include "CelEventEntry.h"
#include "CelFileAndTreeSet.h"
#include "CelFileAndTreeOffset.h"

#include <Rtypes.h>
#include <TString.h>
class RuChain;
class TEventList;

//
// CelEventComponent stores all the information needed to point to a part of an event
// that is located in another data TTree. For a single event this information is
// stored on three TTrees :
// *) the event information: entry index and tree index is stored in the first cel tree
// *) the set information: lists of data tree names is stored in the second cel tree
// *) the set offset: number of entries in previous sets, useful when chaining
//
// The information to navigate between these cel trees is stored by the
// CompositeEventList class, since it is shared by all the components in a skim.
//

class CelEventComponent
 {

  public :

    // construction
    CelEventComponent() ; // needed for ROOT
    CelEventComponent( const TString & componentName ) ;
    virtual ~CelEventComponent() ;

    // Cel Tree Manipulation
    Int_t makeBranches( TTree * entryTree, TTree * fileTree, TTree * offsetTree =0, Int_t bufsize =32000 ) const ;
    Int_t attachToTree( TTree * entryTree, TTree * fileTree, TTree * offsetTree =0 ) ;
  
    // Write
    void registerEntry( TTree & ) ;
    void registerEntry( const TString & filePath, const TString & treeName, Long64_t treeNbEntries, Long64_t entryIndex ) ; // shallow flavor
    void nextSet() ;
  
    /// Shallow read
    Bool_t addSetToChain( RuChain * & chain ) ; 
    Long64_t currentIndexInChain() const ;

    /// Deep read
    void setDataAddress
     ( const TString & branchName,
       void * address ) ;
    void resetSet() ;
    Int_t deepRead() ;

    // USEFUL ?
    TTree * getTree() const ;
    const CelEventEntry & currentEntryIndex() const { return _currentEntryIndex ; }
    const CelFileAndTreeSet & currentFileSet() const { return _currentSet ; }
    CelFileAndTreeOffset & currentFileOffset() { return _currentOffset ; }
    const TString & componentName() const { return _componentName ; }

    // Printing Functions
    // print a single event
    void printEventInfo(const char* options ) const;
    // print all the trees used by this Components
    void printSetInfo(const char* options, const char* prefix ="" ) const;

   private :
  
     // utilities

     // disable copying and assignment
     CelEventComponent( const CelEventComponent & ) ;
     CelEventComponent & operator=( const CelEventComponent & ) ;
  
     // data for shallow read
     TString _componentName ;
     CelEventEntry _currentEntryIndex ;
     CelFileAndTreeSet _currentSet ;
     CelFileAndTreeOffset _currentOffset ;
     
     // additionnal data for deep read
     TString _mainBranchName ;
     void *  _data ;
     TTree * _tree ;

     ClassDef(CelEventComponent,0)

 } ;


#endif

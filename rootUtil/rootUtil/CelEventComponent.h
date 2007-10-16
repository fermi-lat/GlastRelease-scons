
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

#include "CelFileAndTreeSet.h"
#include "CelEventEntry.h"

#include <Rtypes.h>
#include <TString.h>
class TChain;
class TEventList;

//
// CelEventComponent stores all the information needed to point to a part of an event
// that is located in another data TTree. For a single event this information is
// stored on two TTrees :
// *) the event information: entry index and tree index is stored in the first cel tree
// *) the collection information: lists of data tree names is stored in the second cel tree
//
// The information to navigate between these two cel trees is stored by the
// CompositeEventList class, since it is shared by all the components in a skim.
//

class CelEventComponent
 {

  public :

    // construction
    CelEventComponent() ; // needed for ROOT
    CelEventComponent( const TString & componentName ) ;
    virtual ~CelEventComponent() ;

  // Methods and functions
  // set the current event
  void set(TTree& tree);
  // read an event
  Int_t read();
  // Get the Tree that is being read
  TTree* getTree() const ;
  
    // Cel Tree Manipulation
    Int_t makeBranches( TTree & entryTree, TTree & fileTree, Int_t bufsize =32000 ) const ;
    Int_t attachToTree( TTree & entryTree, TTree & fileTree ) ;
  
    /// READING INTERFACE
    // Building a TChain
    Bool_t addToChain( TChain * & chain ) ; 
    Long64_t indexInCurrentSet() const ;

    // USEFUL ? Access
    const CelEventEntry & currentEntryIndex() const { return _currentEntryIndex ; }
    const CelFileAndTreeSet & currentFileSet() const { return _currentSet ; }
    inline const TString & componentName() const { return _componentName ; }

  // Printing Functions
  // print a single event
  void dumpEvent(const char* options) const;
  // print all the trees used by this Components
  void listTrees(const char* options) const;

 

   private :
  
     // disable copying and assignment
     CelEventComponent( const CelEventComponent & ) ;
     CelEventComponent & operator=( const CelEventComponent & ) ;
  
     // data
     TString _componentName ;
     CelEventEntry _currentEntryIndex ;
     CelFileAndTreeSet _currentSet ;

     ClassDef(CelEventComponent,0)

 } ;


#endif

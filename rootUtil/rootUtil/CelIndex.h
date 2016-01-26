
#ifndef CelIndex_h
#define CelIndex_h

/*
* Project: GLAST
* Package: rootUtil
*    File: $Id$
* Authors:
*   EC, Eric Charles, SLAC, echarles@slac.stanford.edu
*   DC, David Chamont, LLR, chamont@llr.in2p3.fr
*
* Copyright (c) 2007
*                   Regents of Stanford University. All rights reserved.
*
*/

class CompositeEventList ;
class CelEventComponent ;

class TTree ;
class TChain ;

#include <TVirtualIndex.h>
#include <RVersion.h> // needed so to have macro ROOT_VERSION_CODE

#include <string>

//
// A way to read the events described in a CEL is :
// 1) make a main TChain from the CEL internal "CompositeEvents" tree
// 2) for each component, make a friend TChain with all corresponding input files added
// 3) in each friend TChain, add an instance oif TVirtualIndex, able to say,
//    for a given entry in the main TChain, what the corresponding entry
//    in the current friend TChain
//
// CelIndex is the subclass of TVirtualIndex which makes the scenario above possible.
// The whole point of this class is to override the method
// TVirtualIndex::GetEntryNumberFriend(), whiwch receive the master TChain as input,
// and must compute the corresponding index in the current friend TChain.
//
// Actually, this index is "PreviousSetsOffset + PreviousTreesOffset + EventIndex"
// * PreviousSetsOffset is the offset from all the previous set of files and trees
// * PreviousTreesOffset is the offset from all the previous TTrees in the current set
// * EventIndex is the Index of the event in the tree that it actually lives on.
// 


class CelIndex : public TVirtualIndex
 {
  public :
  
    // The real useful stuff
    CelIndex( CompositeEventList &, const TString & componentName, TChain * componentChain ) ;
    virtual ~CelIndex() ;
    virtual Long64_t GetEntryNumberFriend( const TTree * celCompositeEvents ) ; // the main method
    virtual Long64_t GetN() const ; // total number of Events in the index

    // Dummy implementation
    CelIndex() ; // Needed for ROOT
#if ROOT_VERSION(5,18,0) <= ROOT_VERSION_CODE
      virtual void UpdateFormulaLeaves( const TTree * ) { return ; }
#else
      virtual void UpdateFormulaLeaves() { return ; }
#endif
    virtual void SetTree( const TTree * tree )
     { fTree = const_cast<TTree*>(tree) ; }
  
    // Unsupported part of the TVirtualIndex interface
    virtual Long64_t GetEntryNumberWithBestIndex( Long64_t, Long64_t ) const
     { MayNotUse("GetEntryNumberWithBestIndex") ; return 0 ; }
    virtual Long64_t GetEntryNumberWithIndex( Long64_t, Long64_t ) const
     { MayNotUse("GetEntryNumberWithIndex") ; return 0; }
    virtual const char * GetMajorName() const
     { MayNotUse("GetMajorName") ; return 0 ; }
    virtual const char * GetMinorName() const
     { MayNotUse("GetMinorName") ; return 0 ; }
    virtual void Append(const TVirtualIndex* /*ind*/ , Bool_t /*delaySort*/) 
      { MayNotUse("Append"); }

  private :
  
    // disable copying and assignment
    CelIndex( const CelIndex & ) ;
    CelIndex & operator=( const CelIndex & ) ;

    // Data
    CompositeEventList * _cel ;      //! Pointer to the CEL
    CelEventComponent * _component ; //! Pointer to the component associated to the current index

//    ClassDef(CelIndex,0)

 } ;


#endif

 

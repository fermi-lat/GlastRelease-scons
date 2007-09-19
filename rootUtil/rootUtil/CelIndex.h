// -*- Mode: c++ -*-
#ifndef CelIndex_h
#define CelIndex_h
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

// Base class headers
#include <TVirtualIndex.h>

// Headers for associated classes

// c++ and stl headers
#include <string>

// forward declares
class TTree;
class CelComponent;
class CompositeEventList;

//
// CelIndex is a sub-class of TVirtualIndex.  It is using a CompositeEventList as an index in to the component trees.
//
// The whole point of this class is to override the method
//   Int_t TVirtualIndex::GetEntryNumberFriend(const TTree* tree);  
//
// Which is called in TTree::LoadTreeFriend(Long64_t entry, TTree* masterTree) to load the friend tree.
// 
// CelIndex::GetEntryNumberFriend(const TTree* tree) goes into the master tree and finds the correct event index
// in the TChain of the component tree.  
//
// The Index in the component tree is
//
//     GlobalOffset + MetaEventOffset + LocalOffset + EventIndex
//
//     GlobalOffset is the offset from starting in the middle of the chain.  Except in special cases it is 0.
//     MetaEventOffset is the offset from all the previous entries in the FileTree
//     LocalOffset is the offset form all the previous TTrees in the current entry in the FileTree
//     EventIndex is the Index of the event in the tree that it actually lives on.
// 


class CelIndex : public TVirtualIndex {

public:
  
  // Build and return a pointer index from a CEL
  static CelIndex* buildIndex(CompositeEventList& skim, const std::string& compontent, TTree* t, Long64_t offset=0);

public:
  
  // ctor's and d'tor
  // Default c'tor.  Needed for ROOT
  CelIndex();
  // Standard c'tor, sets up the index
  CelIndex(CompositeEventList& skim, const std::string& component, TTree* tree, Long64_t offset = 0);
  // D'tor
  virtual ~CelIndex();

  // Methods and functions
  // Get the Entry Number in the CelComponent Tree
  virtual Int_t	GetEntryNumberFriend(const TTree* tree);  
  // Get the total number of Events in the index
  virtual Long64_t GetN() const;

  // Return the offset of the first event
  Long64_t offset() const {
    return _offset;
  }

protected:

  // Pre-Empt or implement stuff from abtract base-class
  virtual void	UpdateFormulaLeaves() { return; }
  virtual void SetTree(const TTree *tree) {
    fTree = const_cast<TTree*>(tree);
  }
  virtual Long64_t	GetEntryNumberWithBestIndex(Int_t /* major */, Int_t /* minor */) const { 
    MayNotUse("GetEntryNumberWithBestIndex"); return 0; 
  };
  virtual Long64_t	GetEntryNumberWithIndex(Int_t /* major */, Int_t /* minor */) const { 
    MayNotUse("GetEntryNumberWithIndex"); return 0; 
  };
  virtual const char*	GetMajorName() const { 
    MayNotUse("GetMajorName"); return 0; 
  };
  virtual const char*	GetMinorName() const { 
    MayNotUse("GetMinorName"); return 0; 
  };

private:
  
  //disable copying and assignment
  CelIndex(const CelIndex& other);
  CelIndex& operator=(const CelIndex& other);

  // Data
  Long64_t     _offset;              //! Global offset (used when Index picks up in the middle of a chain)
  CompositeEventList* _skim;                //! Pointer to skim that does the actual storing
  CelComponent*   _component;           //! Pointer to relevent component of the skim

  ClassDef(CelIndex,0)           // Sub-class of TVirtualIndex implemented using CompositeEventList

};


#endif

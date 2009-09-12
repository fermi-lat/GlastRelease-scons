// -*- Mode: c++ -*-
#ifndef BgDataHandleBase_h
#define BgDataHandleBase_h
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

// Root headers
#include <Rtypes.h>

// forward declares
class TBranch;
class TTree;

//
// BgDataHandleBase is an abstract base class for associating simple data types with 
// root trees.
//
// For example, the sub-classes BgDataHandle<T> manages a datum of type T, and gets it on and off
// TTrees while allowing us to treat the Handle as thought it were of type T itself.
//
// To declare new branches on a TTree
//   makeBranch(TTree& tree, const char* name, Int_t bufsize = 32000) const;
//
// To use a BranchGroup to attach BgDataHandles to branches on an existing TTree
//   attachToTree(TTree& tree, const char* prefix = 0);
//


class BgDataHandleBase {
public:

  // C'tors and D'tor
  // Default c'tor
  BgDataHandleBase(){;}
  // Copy c'tor
  BgDataHandleBase(const BgDataHandleBase& /*other*/){;}
  // D'tor
  virtual ~BgDataHandleBase();

  // Methods & functions
  // Tree Interactions
  virtual TBranch* makeBranch(TTree& tree, const char* name, Int_t bufsize = 32000) const;
  virtual Bool_t attachToTree(TTree& tree, const char* name) const;

  // Type of branch to make,  Returns a Char_t that root uses in TTree::Branch()
  virtual Char_t branchType() const = 0;

protected:

  // The address of the memory where the data is located.  This is passes to Root.
  virtual void* void_ptr() const = 0;

  ClassDef(BgDataHandleBase,0) // Base class for associating simple data types with TTrees
};


#endif

// -*- Mode: c++ -*-
#ifndef TreeRefHandle_h
#define TreeRefHandle_h
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
#include "BranchGroup.h"

// Headers for associated classes
#include "DataHandle.h"

// c++ and stl headers
#include <map>

// forward declares
class TTree;
class TFile;
class TObjArray;
class TArrayL64;

//
// TreeRefHandle stores information needed to point to a part of an event
// that is located in another TTree.
//
// TreeRefHandle keep track of the TTrees where the event components live.
//
// Lists of files and tree names are stored persistently on a TTree by 
// TreeRefHandle.  These are:
//
//  _fileNameList (goes to branch XXX_Files) is the vector of file names
//  _treeNameList (goes to branch XXX_Trees) is the vector of tree names
// 
// At any given time a TreeRefHandle may be keeping track of any number of TTree
//   
// Access to the TTree is by key
//  TTree* getTree(UShort_t key) const ;
//
// The keys are assigned when TTree are first associated with the TreeRefHandle
// 

class TreeRefHandle  : public BranchGroup {

public:

  // c'tors & d'tor
  // Default c'tor.  Needed for ROOT
  TreeRefHandle();
  // Standard c'tor, stores the name of the component 
  TreeRefHandle(const std::string& componentName);
  // D'tor
  virtual ~TreeRefHandle();

  // Methods and functions
  // Reset this object and clear caches
  void reset();
  // Add a tree to the list refered to by this object
  UShort_t addTree(TTree& tree);
  // Get the persistent KEY for a given tree
  UShort_t getKey(TTree* tree) const;
  // Get a given tree using persistent KEY
  TTree* getTree(UShort_t key) const ;
  // Get the Event offset using persistent KEY
  Long64_t getOffset(UShort_t key) const ;

  // Override the methods in BranchGroup.  
  virtual Int_t makeBranches(TTree& tree, const char* prefix = 0, Int_t bufsize = 32000) const;
  virtual Int_t attachToTree(TTree& tree, const char* prefix = 0);

  // Access
  inline UShort_t size() const { return _size; }
  inline Long64_t entries() const { return _entries; }
  inline const std::string& componentName() const { return _componentName; }

  // Printing
  // Print the list of trees
  void show(const char* options = "tf") const;
  // Print information about a single tree
  void printTreeInfo(UShort_t key, const char* options = "tf") const;

protected:

  // Utility function to actually go and get a tree out of a file
  TTree* fetchTree(UShort_t key) const;

private:

  //disable copying and assignment
  TreeRefHandle(const TreeRefHandle& other);
  TreeRefHandle& operator=(const TreeRefHandle& other);

  // data
  std::string            _componentName;   //! Name of the component this object is storing information about
  DataHandle<UShort_t>   _size;            //! Number of trees stored by this TreeRefHandle
  DataHandle<Long64_t>   _entries;         //! Number of events on the trees stored by this TreeRefHandle
  TObjArray*             _treeNameList;    //! Names of trees stored by this TreeRefHandle
  TObjArray*             _fileNameList;    //! Names of files where the TTree live
  TArrayL64*             _treeOffsets;     //! Number of entries in TTrees

  // Cache for fast lookup
  mutable std::map<UShort_t,TTree*> _cache;  //! Maps key -> tree
  mutable std::map<TTree*,UShort_t> _lookup; //! Reverse map, for making sure we don't screw up

  ClassDef(TreeRefHandle,0) // For handling referenece to trees
};


#endif

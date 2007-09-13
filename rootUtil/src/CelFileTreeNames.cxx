// -*- Mode: c++ -*-
#ifndef CelFileTreeNames_cxx
#define CelFileTreeNames_cxx
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

//
// CelFileTreeNames
//
// CelFileTreeNames stores information needed to point to a part of an event
// that is located in another TTree.
//
// CelFileTreeNames keep track of the TTrees where the event components live.
//
// Lists of files and tree names are stored persistently on a TTree by 
// CelFileTreeNames.  These are:
//
//  _fileNameList (goes to branch XXX_Files) is the vector of file names
//  _treeNameList (goes to branch XXX_Trees) is the vector of tree names
// 
// At any given time a CelFileTreeNames may be keeping track of any number of TTree
//   
// Access to the TTree is by key
//  TTree* getTree(UShort_t key) const ;
//
// The keys are assigned when TTree are first associated with the CelFileTreeNames
// 

// This Class's header
#include "rootUtil/CelFileTreeNames.h"

// c++/stl headers
#include <iostream>
#include <assert.h>

// ROOT Headers
#include <TTree.h>
#include <TFile.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TArrayL64.h>

// Other headers from this package
#include "rootUtil/FileUtil.h"
#include "rootUtil/DataHandle.h"
#include "rootUtil/OptUtil.h"


ClassImp(CelFileTreeNames);

CelFileTreeNames::CelFileTreeNames():
  BranchGroup(),
  _componentName(),
  _size(0,*this,"size"),
  _entries(0,*this,"entries"){
  // Default c'tor.  CelComponent Name is not set
  _treeNameList = new TObjArray;
  _fileNameList = new TObjArray;  
  _treeOffsets = new TArrayL64;
}

CelFileTreeNames::CelFileTreeNames(const std::string& componentName):
  BranchGroup(),
  _componentName(componentName),
  _size(0,*this,"size"),
  _entries(0,*this,"entries"){
  // Standard c'tor.  CelComponent Name is set
  _treeNameList = new TObjArray;
  _fileNameList = new TObjArray;  
  _treeOffsets = new TArrayL64;
}

CelFileTreeNames::~CelFileTreeNames(){
  // Delete stuff
  delete _treeNameList;
  delete _fileNameList;
  delete _treeOffsets;
}

void CelFileTreeNames::reset() {
  // Clear the lists
  //
  // Should be called when switch to new input collection
  _size = 0;
  _treeNameList->Clear();
  _fileNameList->Clear();
  _cache.clear();
  _lookup.clear();
}

UShort_t CelFileTreeNames::addTree(TTree& tree) {
  // Add a new tree to the set of trees this object is looking after
  // 
  const char* tName = tree.GetName();
  TFile* f = FileUtil::getFile(tree);
  if ( 0 == f ) {
    // Already warned.  Just return FileUtil::NOKEY
    return FileUtil::NOKEY;
  }
  const char* fName = f->GetName();
  TObjString* tNameSt = new TObjString(tName);
  TObjString* fNameSt = new TObjString(fName);
  _treeNameList->AddLast(tNameSt);
  _fileNameList->AddLast(fNameSt);  
  UShort_t retVal = _size;
  _size = _size + 1;
  // should do something better here  
  _treeOffsets->Set(_size);
  _treeOffsets->AddAt(_entries,retVal);
  _entries = _entries + tree.GetEntries();
  _cache[retVal] = &tree;
  _lookup[&tree] = retVal;
  return retVal;
}

TTree* CelFileTreeNames::getTree(UShort_t key) const {
  // Get a given tree using persistent KEY
  //
  // Return NULL silently if key == FileUtil::NOKEY
  // Warns and returns NULL if tree is not found
  if ( key == FileUtil::NOKEY ) return 0;
  std::map<UShort_t,TTree*>::iterator itrFind = _cache.find(key);
  TTree* tree = itrFind == _cache.end() ? fetchTree(key) : itrFind->second;
  if ( 0 == tree ) {
    // Already warned when fetchTree failed.  More?
    ;
  }
  return tree;
}

UShort_t CelFileTreeNames::getKey(TTree* tree) const {
  // Get the persistent KEY for a given tree
  //
  // Return FileUtil::NOKEY if 'tree' is NULL
  // Warns and returns FileUtil::NOKEY if tree is not found in lookup table
  if ( 0 == tree ) return FileUtil::NOKEY;
  std::map<TTree*,UShort_t>::iterator itrFind = _lookup.find(tree); 
  if ( itrFind == _lookup.end() ) {
    // No found, warn?
    return FileUtil::NOKEY;
  }
  return itrFind->second;
}

Long64_t CelFileTreeNames::getOffset(UShort_t key) const {
  // Get the Event offset using persistent KEY
  //
  // Return 0 if "key" is FileUtil::NOKEY
  // Returns -1 if key is not found
  if ( key == FileUtil::NOKEY ) return 0;
  if ( key >= _size ) return -1;
  return _treeOffsets->At(key);
}

Int_t CelFileTreeNames::makeBranches(TTree& tree, const char* prefix, Int_t bufsize) const {
  // Builds branches on 'tree'
  //
  // The branches will be:
  //   <prefix>size  -> Number of TTree used by this component
  //   <prefix>Tree  -> TObjArray of TObjString with the names of the TTrees
  //   <prefix>File  -> TObjArray of TObjString with the names files where the TTrees live

  // First the base class
  Int_t bVal = BranchGroup::makeBranches(tree,prefix,bufsize);
  if ( bVal < 0 ) return bVal;
  // The list of trees
  std::string tbName; if ( prefix != 0 ) tbName += prefix;
  tbName += "Tree";
  TBranch* bTree = tree.Bronch(tbName.c_str(),_treeNameList->ClassName(),(void*)(&_treeNameList),bufsize,0);
  if ( 0 == bTree ) return -1;
  bVal++;
  // The list of files
  std::string fbName; if ( prefix != 0 ) fbName += prefix;
  fbName += "File";
  TBranch* bFile = tree.Bronch(fbName.c_str(),_fileNameList->ClassName(),(void*)(&_fileNameList),bufsize,0);
  if ( 0 == bFile ) return -1;
  bVal++;
  // The tree sizes
  std::string sbName; if ( prefix != 0 ) sbName += prefix;
  sbName += "Offset";
  TBranch* bSize = tree.Bronch(sbName.c_str(),_treeOffsets->Class_Name(),(void*)(&_treeOffsets),bufsize,0);
  if ( 0 == bSize ) return -1;
  bVal++;
  return bVal;
}

Int_t CelFileTreeNames::attachToTree(TTree& tree, const char* prefix) {
  // Attachs to branches on 'tree'
  //
  // The branches will be:
  //   <prefix>size  -> Number of TTree used by this component
  //   <prefix>Tree  -> TObjArray of TObjString with the names of the TTrees
  //   <prefix>File  -> TObjArray of TObjString with the names files where the TTrees live

  // get rid of old stuff
  reset();
  // First the base class
  Int_t bVal = BranchGroup::attachToTree(tree,prefix);
  if ( bVal < 0 ) return bVal;  
  // The list of trees
  std::string tbName; if ( prefix != 0 ) tbName += prefix;
  tbName += "Tree";
  TBranch* bTree = tree.GetBranch(tbName.c_str());
  if ( 0 == bTree ) return -1;
  bTree->SetAddress((void*)(&_treeNameList));
  bVal++;
  // The list of files
  std::string fbName; if ( prefix != 0 ) fbName += prefix;
  fbName += "File";
  TBranch* bFile = tree.GetBranch(fbName.c_str());
  if ( 0 == bFile ) return -1;
  bFile->SetAddress((void*)(&_fileNameList));  
  bVal++;
  // The tree sizes
  std::string sbName; if ( prefix != 0 ) sbName += prefix;
  sbName += "Offset";
  TBranch* bSize = tree.GetBranch(sbName.c_str());
  if ( 0 == bSize ) return -1;
  bSize->SetAddress((void*)(&_treeOffsets));  
  bVal++;
  return bVal;
}

void CelFileTreeNames::show(const char* options) const {
  // Print the list of trees, one per line
  //
  // if "options" includes 'f' this will print the file names
  // if "options" includes 't' this will print the tree names
  // if "options" includes 'o' this will print the offset in events
  for (UShort_t i(0); i < _size; i++ ) {
    std::cout << "Tree " << i << ":\t";
    printTreeInfo(i,options);    
    std::cout << std::endl;    
  }
}

void CelFileTreeNames::printTreeInfo(UShort_t key, const char* options) const {  
  // Print information about tree with "key" on one line
  //
  // if "options" includes 'f' this will print the file name
  // if "options" includes 't' this will print the tree name
  // if "options" includes 'o' this will print the offset in events
  if (OptUtil::has_option(options,'f')) {
    std::cout <<  _fileNameList->UncheckedAt(key)->GetName() << ' ';
  }
  if (OptUtil::has_option(options,'t')) {
    std::cout << _treeNameList->UncheckedAt(key)->GetName() << ' ';
  }
  if (OptUtil::has_option(options,'o')) {
    std::cout << _treeOffsets->At(key) << ' ';
  }
}

TTree* CelFileTreeNames::fetchTree(UShort_t key) const {
  // Utility function to actually go and get a tree out of a file
  //
  // Return NULL silently if key == FileUtil::NOKEY
  // Warns and returns NULL if tree is not found
  if ( key == FileUtil::NOKEY ) return 0;
  assert( key < _size );
  const char* fName = _fileNameList->UncheckedAt(key)->GetName();

  TFile* f = FileUtil::openFile(fName);
  if ( 0 == f ) { 
    // ROOT has already warned.  Just exit.
    return 0;
  }

  const char* tName = _treeNameList->UncheckedAt(key)->GetName();
  TObject* obj = f->Get(tName);
  
  TTree* tree = dynamic_cast<TTree*>(obj);  
  if ( 0 == tree ) { 
    std::cerr << "Could not find TTree called " << tName << " in file " << fName << std::endl;
    return 0;
  }

  _lookup[tree] = key;
  _cache[key] = tree;
  return tree;
}



#endif

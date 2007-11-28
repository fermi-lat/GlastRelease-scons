
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


#include "rootUtil/CelFileAndTreeSet.h"

// c++/stl headers
#include <iostream>
#include <assert.h>

// ROOT Headers
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TArrayL64.h>

// Other headers from this package
#include "rootUtil/FileUtil.h"
#include "rootUtil/BgDataHandle.h"
#include "rootUtil/OptUtil.h"


ClassImp(CelFileAndTreeSet);

CelFileAndTreeSet::CelFileAndTreeSet()
 : _componentName(),
   _setSize(0,*this,"Set_Size"),
   _treesSize(0,*this,"Set_Entries")
 {
  _treeNames = new TObjArray ;
  _fileNames = new TObjArray ;  
  _treeOffsets = new TArrayL64 ;
 }

CelFileAndTreeSet::CelFileAndTreeSet(const TString & componentName )
 : _componentName(componentName),
   _setSize(0,*this,"Set_Size"),
   _treesSize(0,*this,"Set_Entries")
 {
  _fileNames = new TObjArray ;  
  _treeNames = new TObjArray ;
  _treeOffsets = new TArrayL64 ;
 }

CelFileAndTreeSet::~CelFileAndTreeSet()
 {
  delete _fileNames ;
  delete _treeNames ;
  delete _treeOffsets ;
 }

void CelFileAndTreeSet::reset()
 {
  // Clear the lists
  // Should be called when switch to new input collection
  _setSize = 0 ;
  _fileNames->Clear() ;
  _treeNames->Clear() ;
  _cache.clear() ;
  _lookup.clear() ;
}

UShort_t CelFileAndTreeSet::addTree(TTree& tree) {
  // Add a new tree to the set of trees this object is looking after
  // 
  const char* tName = tree.GetName();
  TFile* f = FileUtil::getFile(tree);
  if ( 0 == f ) {
    // Already warned.  Just return FileUtil::NOKEY
    return FileUtil::NOKEY;
  }
  const char* fName = f->GetName();
  TObjString* fNameSt = new TObjString(fName);
  TObjString* tNameSt = new TObjString(tName);
  _fileNames->AddLast(fNameSt);  
  _treeNames->AddLast(tNameSt);
  UShort_t retVal = _setSize;
  _setSize = _setSize + 1;
  // should do something better here  
  _treeOffsets->Set(_setSize);
  _treeOffsets->AddAt(_treesSize,retVal);
  _treesSize = _treesSize + tree.GetEntries();
  _cache[retVal] = &tree;
  _lookup[&tree] = retVal;
  return retVal;
}

TTree* CelFileAndTreeSet::getTree(UShort_t key) const {
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

UShort_t CelFileAndTreeSet::getKey( TTree * tree ) const
 {
  // Get the persistent KEY for a given tree
  // Return FileUtil::NOKEY if 'tree' is NULL
  // Warns and returns FileUtil::NOKEY if tree is not found in lookup table
  if ( 0 == tree ) return FileUtil::NOKEY;
  std::map<TTree*,UShort_t>::iterator itrFind = _lookup.find(tree) ; 
  if ( itrFind == _lookup.end() )
   {
    // No found, warn?
    return FileUtil::NOKEY ;
   }
  return itrFind->second ;
 }

Long64_t CelFileAndTreeSet::getOffset(UShort_t key) const {
  // Get the Event offset using persistent KEY
  //
  // Return 0 if "key" is FileUtil::NOKEY
  // Returns -1 if key is not found
  if ( key == FileUtil::NOKEY ) return 0;
  if ( key >= _setSize ) return -1;
  return _treeOffsets->At(key);
}


//======================================================================================
//
// BranchGroup methods
//
//======================================================================================


Int_t CelFileAndTreeSet::makeBranches( TTree & tree, const char * prefix, Int_t bufsize ) const
 {
  // First the base class
  Int_t bVal = BranchGroup::makeBranches(tree,prefix,bufsize) ;
  if ( bVal < 0 ) return bVal ;
  
  std::string prefixName = "" ;
  if ( prefix != 0 )
   { prefixName = prefix ; prefixName  += "_" ; }
  
  // The list of trees
  std::string tbName = prefixName + "Set_TreeNames";
  TBranch* bTree = tree.Bronch(tbName.c_str(),_treeNames->ClassName(),(void*)(&_treeNames),bufsize,0);
  if ( 0 == bTree ) return -1;
  bVal++;
  
  // The list of files
  std::string fbName = prefixName + "Set_FileNames";
  TBranch* bFile = tree.Bronch(fbName.c_str(),_fileNames->ClassName(),(void*)(&_fileNames),bufsize,0);
  if ( 0 == bFile ) return -1;
  bVal++;
  
  // The tree sizes
  std::string sbName = prefixName + "Set_TreeOffsets";
  TBranch* bSize = tree.Bronch(sbName.c_str(),_treeOffsets->Class_Name(),(void*)(&_treeOffsets),bufsize,0);
  if ( 0 == bSize ) return -1;
  bVal++;
  return bVal;
}

Int_t CelFileAndTreeSet::attachToTree( TTree & tree, const char * prefix )
 {
  reset() ;
  Int_t bVal = BranchGroup::attachToTree(tree,prefix) ;
  if ( bVal < 0 ) return bVal ;  
  std::string prefixName = "" ;
  if ( prefix != 0 )
   { prefixName = prefix ; prefixName  += "_" ; }
  
  // The list of trees
  std::string tbName = prefixName + "Set_TreeNames" ;
  TBranch* bTree = tree.GetBranch(tbName.c_str()) ;
  if ( 0 == bTree ) return -1 ;
  bTree->SetAddress((void*)(&_treeNames)) ;
  bVal++ ;
  
  // The list of files
  std::string fbName = prefixName + "Set_FileNames";
  TBranch* bFile = tree.GetBranch(fbName.c_str());
  if ( 0 == bFile ) return -1 ;
  bFile->SetAddress((void*)(&_fileNames)) ;  
  bVal++ ;
  
  // The tree offsets
  std::string obName = prefixName + "Set_TreeOffsets";
  TBranch* bSize = tree.GetBranch(obName.c_str()) ;
  if ( 0 == bSize ) return -1 ;
  bSize->SetAddress((void*)(&_treeOffsets)) ;  
  bVal++ ;
  
  return bVal ;
 }



//======================================================================================
//
// Reading Interface
//
//======================================================================================


Bool_t CelFileAndTreeSet::addToChain( TChain * & chain )
 {
  UShort_t treeIndex ;
  for ( treeIndex=0 ; treeIndex<_setSize ; treeIndex++ )
   {
	const char * tName = _treeNames->UncheckedAt(treeIndex)->GetName() ;
	if (chain==0)
     { chain = new TChain(tName,"Cel Component Chain") ; }
    const char * fName = _fileNames->UncheckedAt(treeIndex)->GetName() ;
    chain->AddFile(fName,-1,tName) ;
   }
  return kTRUE ;
 }



//======================================================================================
//
// Printing & debugging
//
//======================================================================================


// Print the list of trees, one per line
void CelFileAndTreeSet::printTreesInfo( const char * options, const char * prefix ) const
 {
  UShort_t i ;
  for ( i=0 ; i<_setSize; i++ )
   {
    std::cout << "Tree " << i << ":\t" ;
    printTreeInfo(i,options) ;    
    std::cout << std::endl ;    
   }
 }

// Print information about tree i
// if "options" includes 'f' this will print the file name
// if "options" includes 't' this will print the tree name
// if "options" includes 'o' this will print the offset in events
void CelFileAndTreeSet::printTreeInfo( UShort_t treeIndex, const char * options, const char * prefix ) const
 {  
  if (OptUtil::has_option(options,'f'))
   { std::cout <<  _fileNames->UncheckedAt(treeIndex)->GetName() << ' '; }
  if (OptUtil::has_option(options,'t'))
   { std::cout << _treeNames->UncheckedAt(treeIndex)->GetName() << ' '; }
  // why not "Unchecked" below ??
  if (OptUtil::has_option(options,'o'))
   { std::cout << _treeOffsets->At(treeIndex) << ' '; }
 }



//======================================================================================
//
// Utilities
//
//======================================================================================


// Utility function to actually go and get a tree out of a file
// Return 0 silently if treeIndex == FileUtil::NOKEY
// Warns and returns 0 if tree is not found
TTree * CelFileAndTreeSet::fetchTree( UShort_t treeIndex ) const
 {
  if ( treeIndex == FileUtil::NOKEY ) return 0 ;
  assert( treeIndex < _setSize ) ;
  
  // open file
  const char * fName = _fileNames->UncheckedAt(treeIndex)->GetName() ;
  TFile * f = FileUtil::openFile(fName) ;
  if ( 0 == f )
   { 
    // ROOT has already warned.  Just exit.
    return 0 ;
   }

  // get tree
  const char * tName = _treeNames->UncheckedAt(treeIndex)->GetName() ;
  TObject * obj = f->Get(tName) ;
  TTree * tree = dynamic_cast<TTree*>(obj) ;  
  if ( 0 == tree )
   { 
    std::cerr
      << "[CelFileAndTreeSet::fetchTree] "
      << "Could not find TTree called " << tName
      << " in file " << fName
      << std::endl ;
    return 0 ;
   }

  _lookup[tree] = treeIndex ;
  _cache[treeIndex] = tree ;
  return tree ;
 }




// -*- Mode: c++ -*-
#ifndef BranchGroup_cxx
#define BranchGroup_cxx
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

// BranchGroup
//
// BranchGroup is a abstract class to manage a group set of values that
// get written to a TTree
//
// Any class that uses BgDataHandle to interact with TTree should inherit 
// from BranchGroup
//
// To use a BranchGroup to declare new branches on a TTree
//   makeBranches(TTree& tree, const char* prefix = 0, Int_t bufsize = 32000) 
//
// To use a BranchGroup to attach BgDataHandles to branches on an existing TTree
//   attachToTree(TTree& tree, const char* prefix = 0);
//
// In either case the branches will have names "<prefix><branchName>"
//

// This Class's header
#include "rootUtil/BranchGroup.h"

// c++/stl headers
#include <iostream>
#include <assert.h>

// Other headers from this package
#include "rootUtil/BgDataHandle.h"

ClassImp(BranchGroup);


BranchGroup::~BranchGroup(){
  // D'tor is a no-op
  ;
}


// Builds branches with names 'prefix_'branchName' on 'tree'
// returns the number of newly build branches, or -1 for failure
Int_t BranchGroup::makeBranches( TTree & tree, const char * prefix, Int_t bufsize ) const
 {
  std::string prefixName ;
  if ( prefix != 0 )
   { prefixName = prefix ; prefixName  += "_" ; }
  
  Int_t nb(0) ;
  std::string fullName ;
  std::list<std::string>::const_iterator itr ;
  for ( itr = _branchNameList.begin() ;	itr != _branchNameList.end() ; itr++ )
   {
    BgDataHandleBase* branch = _branchMap.find(*itr)->second;
    assert( 0 != branch ) ;
    fullName = prefixName + *itr ;
    if ( branch->makeBranch(tree,fullName.c_str(),bufsize) )
     { nb += 1 ; }
    else
     { return -1 ; }
   }
  return nb ;
 }

// Attaches to branches with names 'prefix''branchName' on 'tree'
// returns the number of newly attached branches, or -1 for failure
Int_t BranchGroup::attachToTree( TTree & tree, const char * prefix ) const
 {
  std::string prefixName ;
  if ( prefix != 0 )
   { prefixName = prefix ; prefixName  += "_" ; }

  Int_t nb(0) ;
  std::string fullName ;
  std::list<std::string>::const_iterator itr ;
  for ( itr = _branchNameList.begin() ;	itr != _branchNameList.end() ; itr++ )
   {
    BgDataHandleBase* branch = _branchMap.find(*itr)->second;
    assert( 0 != branch );
    fullName = prefixName + *itr ;
    if ( branch->attachToTree(tree,fullName.c_str()) )
     { nb += 1; }
    else
     { return -1 ;  }
  }
  return nb ; 
 }

Bool_t BranchGroup::addBranch(const char* name, BgDataHandleBase& branch) {
  // Add a branch to this group.  The branch will be stored under the 'name'
  std::string theName(name);
  std::map<std::string,BgDataHandleBase*>::iterator itrFind = _branchMap.find(theName);
  if ( itrFind != _branchMap.end() ) {
    return kFALSE;
  }
  _branchNameList.push_back(theName);
  _branchMap[theName] = &branch;
  return kTRUE;
}


BgDataHandleBase* BranchGroup::getBranch(const char* name) {
  // Get a branch by name
  std::string theName(name);
  std::map<std::string,BgDataHandleBase*>::iterator itrFind = _branchMap.find(theName);
  return itrFind == _branchMap.end() ? 0 : itrFind->second;
}

#endif

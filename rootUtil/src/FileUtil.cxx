// -*- Mode: c++ -*-
#ifndef FileUtil_cxx
#define FileUtil_cxx
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

// This Class's header
#include "rootUtil/FileUtil.h"

// c++/stl headers
#include <string>
#include <iostream>

// ROOT Headers
#include <TROOT.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>

  
TFile* FileUtil::openFile(const char* fileName) {
  // Open a file.  Return NULL if failed
  TFile* f = TFile::Open(fileName);
  if ( 0 == f ) {
    // ROOT already warned.  Just return NULL 
    return 0;
  }
  return f;  
}


TFile* FileUtil::getFile(TTree& tree) {
  // Get the TFile associated to a TTree.  Return NULL if failed
  //
  // This can fail if the TTree has not yet been associated with a TFile.
  // For example if you are making a composite event list at the same tile as the 
  // files it points to.
  //
  // There are two ways to solve that:
  //   Use TFile::cd() to switch to a directory before making the TTrees you point to.
  //   Use TTree::SetDirectory(aTFile) to set the TTree's directory
  //
  // Either case implies that you have to select and open the output file with the 
  // event component before trying to make a composite event list.
  TFile* f = tree.GetCurrentFile();
  if ( 0 == f ) {
    // Tree doesn't have an associated file yet.  
    // WARN and return NULL
    std::cerr << "TTree called " << tree.GetName() << " at " << &tree 
	      << " does not have an associated file." << std::endl
	      << "You have to associated TTree with TFiles before you can persist "  << std::endl
	      << "pointers to the event components they contain" << std::endl
	      << "If you are making the CompositeEventList in the same job as the TTree it " << std::endl
	      << "points to you need to do one of two things" << std::endl
	      << "\t1)Use TFile::cd() to switch to a directory before making the TTree" << std::endl
	      << "\t2)Use TTree::SetDirectory(aTFile) to set the TTree's directory" << std::endl;
    return 0;
  }
  // This is a TFile of some kind.  Good enough.
  return f;
}  

#endif

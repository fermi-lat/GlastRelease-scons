
/*
* Project: GLAST
* Package: rootUtil
* File: RuUtil.cxx
* Authors:
*   DC, David Chamont, LLR,  chamont@llr.in2p3.fr
*   EC, Eric Charles,  SLAC, echarles@slac.stanford.edu
*/

// Main header
#include "rootUtil/RuUtil.h"

// ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>

// c++/stl headers
#include <iostream>
#include <string>



//====================================================================
// for files
//====================================================================


rootUtil::AutoCd::AutoCd( TDirectory * newDir )
 : _oldDir(gDirectory)
 { if (newDir!=0) newDir->cd() ; }

rootUtil::AutoCd::~AutoCd()
 { if (_oldDir!=0) _oldDir->cd() ; }

TFile * rootUtil::openFile( const char * fileName )
 {
  TFile * f = TFile::Open(fileName) ;
  if ( f == 0 )
   {
    return 0 ;
   }
  return f ;  
 }

TFile * rootUtil::getFile( TTree & tree)
 {
  // Get the TFile associated to a TTree.  Return NULL if failed
  //
  // This can fail if the TTree has not yet been associated with a TFile.
  // For example if you are making a composite event list at the same tile as the 
  // files it points to.
  //
  // There are two ways to solve that:
  //   Use TFile::cd() to switch to a directory before making the TTrees you point to.
  //   Use TTree::SetDirectory(aTFile) to set the TTree's directory
  
  TFile * f = tree.GetCurrentFile() ;
  if ( f == 0 )
   {
    // Tree doesn't have an associated file yet.  
    // WARN and return NULL
    std::cerr
      << " TTree called " << tree.GetName() << " at " << &tree 
      << " does not have an associated file." << std::endl ;
    return 0 ;
   }
  
  // This is a TFile of some kind.  Good enough.
  return f ;
 }  



//====================================================================
// for options
//====================================================================


Bool_t rootUtil::has_option( const char * options, char opt )
 {
  std::string os(options) ;
  return (os.find(opt)!=os.npos?kTRUE:kFALSE) ;
 }  


// -*- Mode: c++ -*-
//#ifndef BgDataHandle_cxx
//#define BgDataHandle_cxx
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

// BgDataHandle<T>
//
// BgDataHandle<T> is a template for associated built-in types with 
// branches on root TTrees.  It is mainly a wrapper around an datum of a built-in type which
// function like the underlying datum
//
// Specifically, BgDataHandle<T> implements the casting operator and assignement operators
//  BgDataHandle<T>::operator T()
//  BgDataHandle<T>
//
// Any class that uses BgDataHandle to interact with TTree should inherit 
// from BranchGroup
//  
// In the BranchGroup sub-class you build a BgDataHandle by with a reference to the BranchGroup.  
// For example:
//   BgDataHandle<Int_t> x;    in the header file
//   x(0,*this,"xvar";       in the source file
// Will cause an Int_t branch called "xvar" to be associated with any TTree that BranchGroup 
// is linked with.
//
//

#include <rootUtil/BgDataHandle.h>
#include <rootUtil/BranchGroup.h>

//ClassImpT(BgDataHandle,T) ;

template <class T>
BgDataHandle<T>::~BgDataHandle()
 {}

// Partial specialiaztion for return branchType
template<> Char_t BgDataHandle<char *>::branchType() const { return 'C' ; }
template<> Char_t BgDataHandle<Bool_t>::branchType() const { return 'O' ; }
template<> Char_t BgDataHandle<Char_t>::branchType() const { return 'B' ; } 
template<> Char_t BgDataHandle<UChar_t>::branchType() const { return 'b' ; } 
template<> Char_t BgDataHandle<Short_t>::branchType() const { return 'S' ; } 
template<> Char_t BgDataHandle<UShort_t>::branchType() const { return 's' ; }
template<> Char_t BgDataHandle<Int_t>::branchType() const { return 'I' ; }
template<> Char_t BgDataHandle<UInt_t>::branchType() const { return 'i' ; }
template<> Char_t BgDataHandle<Long64_t>::branchType() const { return 'L' ; }
template<> Char_t BgDataHandle<ULong64_t>::branchType() const { return 'l' ; } 
template<> Char_t BgDataHandle<Float_t>::branchType() const { return 'F' ; } 
template<> Char_t BgDataHandle<Double_t>::branchType() const { return 'D' ; }

void BgDataHandleInstance::init() {
  // This is here to force template loading
  static Bool_t done=kFALSE;
  if ( done ) return;
  done = kTRUE;
  static BgDataHandle<char *> ins_string;
  static BgDataHandle<Bool_t> ins_bool;
  static BgDataHandle<Char_t> ins_char;
  static BgDataHandle<UChar_t> ins_uchar;
  static BgDataHandle<Short_t> ins_short;
  static BgDataHandle<UShort_t> ins_ushort;
  static BgDataHandle<Int_t> ins_int;
  static BgDataHandle<UInt_t> ins_uint;
  static BgDataHandle<Long64_t> ins_long;
  static BgDataHandle<ULong64_t> ins_ulong;
  static BgDataHandle<Float_t> ins_float;
  static BgDataHandle<Double_t> ins_double;
  return;
}

//#endif

// -*- Mode: c++ -*-
//#ifndef DataHandle_cxx
//#define DataHandle_cxx
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

// DataHandle<T>
//
// DataHandle<T> is a template for associated built-in types with 
// branches on root TTrees.  It is mainly a wrapper around an datum of a built-in type which
// function like the underlying datum
//
// Specifically, DataHandle<T> implements the casting operator and assignement operators
//  DataHandle<T>::operator T()
//  DataHandle<T>
//
// Any class that uses DataHandle to interact with TTree should inherit 
// from BranchGroup
//  
// In the BranchGroup sub-class you build a DataHandle by with a reference to the BranchGroup.  
// For example:
//   DataHandle<Int_t> x;    in the header file
//   x(0,*this,"xvar";       in the source file
// Will cause an Int_t branch called "xvar" to be associated with any TTree that BranchGroup 
// is linked with.
//
//

// This Class's header
#include "rootUtil/DataHandle.h"

// Other headers in this package
#include "rootUtil/BranchGroup.h"


//ClassImpT(DataHandle,T) ;


// D'tor
template <class T>
DataHandle<T>::~DataHandle(){
  // D'tor is a no-op
  ;
}

// Partial specialiaztion for return branchType
Char_t DataHandle<Bool_t>::branchType() const { 
  // For Bool_t use '0'
  return 'O';
}
Char_t DataHandle<Char_t>::branchType() const { 
  // For Char_t use 'B'
  return 'B';
} 
Char_t DataHandle<UChar_t>::branchType() const { 
  // For UChar_t use 'b'
  return 'b';
} 
Char_t DataHandle<Short_t>::branchType() const { 
  // For Short_t use 'S'
  return 'S';
} 
Char_t DataHandle<UShort_t>::branchType() const { 
  // For UShort_t use 's'
  return 's';
}
Char_t DataHandle<Int_t>::branchType() const { 
  // For Int_t use 'I'
  return 'I';
}
Char_t DataHandle<UInt_t>::branchType() const { 
  // For UInt_t use 'i'
  return 'i';
}
Char_t DataHandle<Long64_t>::branchType() const { 
  // For Long64_t use 'L'
  return 'L';
}
Char_t DataHandle<ULong64_t>::branchType() const { 
  // For ULong64_t use 'l'
  return 'l';
} 
Char_t DataHandle<Float_t>::branchType() const { 
  // For Float_t use 'F'
  return 'F';
} 
Char_t DataHandle<Double_t>::branchType() const { 
  // For Double_t use 'D'
  return 'D';
}

void DataHandleInstance::init() {
  // This is here to force template loading
  static Bool_t done=kFALSE;
  if ( done ) return;
  done = kTRUE;
  static DataHandle<Bool_t> ins_bool;
  static DataHandle<Char_t> ins_char;
  static DataHandle<UChar_t> ins_uchar;
  static DataHandle<Short_t> ins_short;
  static DataHandle<UShort_t> ins_ushort;
  static DataHandle<Int_t> ins_int;
  static DataHandle<UInt_t> ins_uint;
  static DataHandle<Long64_t> ins_long;
  static DataHandle<ULong64_t> ins_ulong;
  static DataHandle<Float_t> ins_float;
  static DataHandle<Double_t> ins_double;
  return;
}

//#endif

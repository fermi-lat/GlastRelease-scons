// -*- Mode: c++ -*-
#ifndef DataHandle_h
#define DataHandle_h
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
#include "DataHandleBase.h"

#include <assert.h>

// Headers for associated classes
#include "BranchGroup.h"

// forward declares
class TBranch;
class TTree;


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


template <class T>
class DataHandle : public DataHandleBase {

public:
  
  // c'tors & d'tor
  // Default c'tor.  Needed for ROOT
  DataHandle():DataHandleBase(){;}
  // Standard c'tor.  Used to add this to a BranchGroup
  DataHandle(const T& val, BranchGroup& group, const char* name):DataHandleBase(),_value(val){
    Bool_t ok = group.addBranch(name,*this);
    assert(ok);
  }
  // D'tor
  virtual ~DataHandle();
  
  // Operators
  // Assignment from T
  DataHandle<T>& operator=(const T& other) { _value = other; return *this; }
  
  // Cast to T
  //operator const T() const { return _value; }
  operator T() const { return _value; }
  
  // Access and setting
  // Type of branch to make.  Returns a Char_t that root uses in TTree::Branch()
  virtual Char_t branchType() const;

protected:

  // Return address of data
  const T* const_ptr() const { return &_value; }
  T* ptr() { return &_value; }
  virtual void* void_ptr() const { return (void*)&_value; }

private:

  //disable copying
  DataHandle(const DataHandle<T>& other);

  // Data 
  T _value;                  // The Datum being managed

  ClassDefT(DataHandle<T>,0) // Type-safe template for handling simple types
};

// Even more special ROOT declarations for template classes:
ClassDefT2(DataHandle,T) ;

#ifndef __CINT__
#ifdef GLAST_COMP_INST
#include "rootUtil/../DataHandle.cc"
#endif // GLAST_COMP_INST
#endif // __CINT__


//
// This is here to force template loading
namespace DataHandleInstance {
  void init();
}

#endif

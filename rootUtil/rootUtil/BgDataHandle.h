// -*- Mode: c++ -*-
#ifndef BgDataHandle_h
#define BgDataHandle_h
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
#include "BgDataHandleBase.h"

#include <assert.h>

// Headers for associated classes
#include "BranchGroup.h"

// forward declares
class TBranch;
class TTree;


//
// BgDataHandle<T> is a template for associated built-in types with 
// branches on root TTrees.  It is mainly a wrapper around an datum of a built-in type which
// function like the underlying datum
//
// Specifically, BgDataHandle<T> implements the casting operator and assignement operators
//  BgDataHandle<T>::operator T()
//  BgDataHandle<T>::operator=( const T & )
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


template <class T>
class BgDataHandle : public BgDataHandleBase {

public:
  
  // c'tors & d'tor
  BgDataHandle() : BgDataHandleBase() {;} // Needed for ROOT
  // Standard c'tor.  Used to add this to a BranchGroup
  BgDataHandle( const T & val, BranchGroup & group, const char * name )
   : BgDataHandleBase(),_value(val)
   {
    Bool_t ok = group.addBranch(name,*this) ;
    assert(ok);
   }
  virtual ~BgDataHandle() ;
  
  // Operators
  BgDataHandle<T> & operator=( const T & other)
   { _value = other ; return *this ; }
  operator T() const { return _value ; }
  
  // Access and setting
  // Type of branch to make.  Returns a Char_t that root uses in TTree::Branch()
  virtual Char_t branchType() const ;

protected:

  // Return address of data
  const T * const_ptr() const { return &_value ; }
  T * ptr() { return &_value ; }
  virtual void * void_ptr() const { return (void*)&_value ; }

private:

  // disable copy construction
  BgDataHandle( const BgDataHandle<T> & other) ;

  // Data 
  T _value ;                  // The Datum being managed

  ClassDefT(BgDataHandle<T>,0) // Type-safe template for handling simple types
} ;

// Even more special ROOT declarations for template classes:
ClassDefT2(BgDataHandle,T) ;

#ifndef __CINT__
#ifdef GLAST_COMP_INST
#include "rootUtil/../BgDataHandle.cc"
#endif // GLAST_COMP_INST
#endif // __CINT__


//
// This is here to force template loading
namespace BgDataHandleInstance {
  void init();
}

#endif

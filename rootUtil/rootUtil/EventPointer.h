// -*- Mode: c++ -*-
#ifndef EventPointer_h
#define EventPointer_h
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

// forward declares
class TreeRefHandle;

//
// EventPointer stores information needed to point to a part of an event
// that is located in another TTree.
//
// EventPointer uses TreeRefHandle to associated TTrees with UShort_t keys.
//
// EventPointer stores two pieces of information for each event.
//
//  _eventIndex (goes to branch XXX_EvtIndex)  is the index of the event component in the TTree
//  _treeKey    (goes to branch XXX_TreeIndex) is the index of the TTree in the associated TreeRefHandle
//


class EventPointer :public BranchGroup {

public:

  // c'tors & d'tor
  // Default c'tor.  Needed for ROOT
  EventPointer(); 
  // Standard c'tor, stores the name of the component 
  EventPointer(const std::string& componentName);
  // Copy c'tor
  EventPointer(const EventPointer& other);
  // D'tor
  virtual ~EventPointer();

  // Operators
  // Assignment operator
  EventPointer& operator=(const EventPointer& other);

  // Methods and functions
  // set the current event
  void set(TTree& tree, TreeRefHandle& handle);
  // read an event
  Int_t read(const TreeRefHandle& handle);
  // get the tree that is being read
  TTree* getTree(const TreeRefHandle& handle) const;

  // Access
  inline const std::string& componentName() const { return _componentName; }
  inline Long64_t eventIndex() const { return _eventIndex; }
  inline UShort_t treeKey() const { return _treeKey; }

  // Printing 
  // print the info about the current event
  void printEventInfo(const char* options) const;  

protected:


private:
  
  // Data
  std::string            _componentName;  //! Name of the component this object is storing information about
  DataHandle<Long64_t>   _eventIndex;     //! Index of the current event in the TTree
  DataHandle<UShort_t>   _treeKey;        //! Key for the current TTree in the associated TreeRefHandle
  
  ClassDef(EventPointer,0)                // Handles Event Level Information to point to Event Components

};


#endif

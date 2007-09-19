// -*- Mode: c++ -*-
#ifndef CelComponent_h
#define CelComponent_h
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

// Root headers
#include <Rtypes.h>

// c++ and stl headers
#include <string>

// Headers of Data Members
#include "CelFileAndTreeNames.h"
#include "CelEntryIndex.h"

// Forward Declares
class TChain;
class TEventList;

//
// CelComponent stores all the information needed to point to a part of an event
// that is located in another TTree.
//
// This for a single event this information is stored on two TTrees. 
//   The event information: index and tree key is stored on one tree
//   The collection information: lists of trees is stored on the other tree
//
// The information to navigate between these two tree is stored by the CompositeEventList class
// since it is shared by all the components in a skim.
//

class CelComponent {

public:

  // ctor's and d'tor
  // Default c'tor.  Needed for ROOT
  CelComponent(); 
  // Standard c'tor, stores the name of the component 
  CelComponent(const std::string& componentName);
  // D'tor
  virtual ~CelComponent();

  // Methods and functions
  // set the current event
  void set(TTree& tree);
  // read an event
  Int_t read();
  // Get the Tree that is being read
  TTree* getTree() const;
  // Tree Manipulation
  Int_t makeBranches(TTree& eventTree, TTree& fileTree, Int_t bufsize = 32000) const;
  Int_t attachToTree(TTree& eventTree, TTree& fileTree);
  // Building a TChain
  Bool_t addToChain(TChain*& chain);  
  // Get the offset to the first event
  Long64_t getLocalOffset() const;
  // Get the index of the event the local chain (ie, in this meta event)
  Long64_t getIndexInLocalChain() const;

  // Access
  const CelEntryIndex& eventPointer() const { return _event; }
  const CelFileAndTreeNames& treePointer() const { return _tree; }
  inline const std::string& componentName() const { return _componentName; }

  // Printing Functions
  // print a single event
  void dumpEvent(const char* options) const;
  // print all the trees used by this Components
  void listTrees(const char* options) const;

 

protected:


private:
  
  //disable copying and assignment
  CelComponent(const CelComponent& other);
  CelComponent& operator=(const CelComponent& other);
  
  // Data
  std::string   _componentName;    //! Name of this component
  CelEntryIndex  _event;            //! Interface with the per-event data
  CelFileAndTreeNames _tree;             //! Interface with the file reference data

  ClassDef(CelComponent,0)            // Class to handle references to a single Event component

};


#endif

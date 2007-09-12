// -*- Mode: c++ -*-
#ifndef CompositeEventList_h
#define CompositeEventList_h
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
#include <TObject.h>

// Headers for associated classes
#include "DataHandle.h"

// c++ and stl headers
#include <map>
#include <vector>

// forward declares
class TTree;
class TFile;
class TChain;
class TVirtualIndex;
class TCollection;
class CelComponent;
class TObjArray;

//
// CompositeEventList associates infomation about events that is stored in several trees.
// 
// Unlike using friend trees in root, this can include sparse collections without 
// requiring you to build a lookup index.
// 
// The information is stored on 3 trees
// 
//   Events:  Tree with 2 Branches per input component, one entry per event
//      Comp_EvtIndex/L  -> Index of Event in the Tree it lives on
//      Comp_TreeIndex/s -> Index of Tree in the Associated vector of trees
//   
//   Files:   Tree with 5 Branches per input component, one entry per input collection
//      Comp_size/s      -> Number of Trees referenced in current entry in File Tree.
//      Comp_entries/L   -> Number of Entries in all trees references in current entry in File Tree
//      Comp_Tree        -> TClonesArray<TObjString> with the tree names (should all be the same)
//      Comp_File        -> TClonesArray<TObjString> with the file names
//      Comp_Offset      -> TArrayL64 with the offsets for each of the Trees in the current entry in File Tree
//
//   Links:  Tree with 3 Branches total, one entry per event
//      Link_EventIndex  -> Absolute index of this event
//      Link_MetaIndex   -> Index of the related entry in the File Tree
//      Link_MetaOffset  -> Offset in events for this entry in the File Tree
// 

class CompositeEventList : public BranchGroup, public TObject {

protected:

  typedef std::pair<TTree*,CelComponent*> TreeAndComponent;

public:

  // ctor's and d'tor
  // Default c'tor.
  CompositeEventList(); 
  // Build a pointer skim from the three trees
  CompositeEventList(TTree& eventTree, TTree& linkTree, TTree& fileTree);
  // D'tor
  virtual ~CompositeEventList();

  // Methods and functions
  // Make a new file.  Will also delare TTree for storing the pointers
  TFile* makeFile(const Char_t* fileName, const Char_t* options);
  // Open an existing file which contains a pointer skim
  TFile* openFile(const Char_t* fileName);

  // Add a component by name.  This is only need when writing.  On read these are discovered.
  UInt_t addComponent(const std::string& name);
  // Set up an event.  Grab the status of a set of TTrees
  Long64_t fillEvent(TObjArray& trees);
  Long64_t fillEvent(std::vector<TTree*>& trees);
  // Store up the List of TTrees that have been used so far
  Long64_t fillMeta();
  // read an event
  Int_t read(Long64_t iEvt);
  // read only the event tree data
  Int_t readEventAddress(Long64_t iEvt);
  // Build the TChain for a component. 
  TChain* buildChain(UInt_t index);
  // Build the TEventList for a component. 
  TVirtualIndex* buildEventIndex(UInt_t index, Long64_t& offset, TTree* tree);
  // Build all the Chains
  TChain* buildLinks(TObjArray* chainList = 0, Bool_t setFriends = kTRUE);


  // Access
  // Number of entries in this skim
  Long64_t entries() const;
  // Get the index of the current entry in the event tree
  Long64_t eventIndex() const { return _eventIndex; }
  // Get the index of the current entry in the file tree
  Long64_t metaIndex() const { return _metaIndex; }
  // Get the offset of the current entry in the file tree
  Long64_t metaOffset() const { return _metaOffset; }
  // Return the tree with the Event component entries
  TTree* eventTree() { return _eventTree; }
  // Return the tree with the links between events and meta data
  TTree* linkTree() { return _linkTree; }
  // Return the tree with the meta data
  TTree* fileTree() { return _fileTree; }
  // Get the number of components
  UInt_t nComponent() const { return _compList.size(); }
  // Get the index of a component by name.  Returns 0xFFFFFFFF if nonesuch component exists
  UInt_t componentIndex(const std::string& name) {
    std::map<std::string,UInt_t>::const_iterator itrFind = _compMap.find(name);
    return itrFind != _compMap.end() ? itrFind->second : 0xFFFFFFFF;
  }
  // Get a Tree that is being read
  TTree* getTree(UInt_t index) const {
    return (index < _compList.size()) ? _compList[index].first : 0;
  }
  // Get an event Compnent that is being read
  CelComponent* getComponent(UInt_t index) const{
    return (index < _compList.size()) ? _compList[index].second : 0;
  }

  // Printing
  // Dump a set of event component pointers and the list of TTree where they live
  void printout(const char* options, UInt_t nEvent=10, UInt_t startEvent=0);
  // Dump a single event
  void dumpEvent(const char* options);
  // Dump the list of TTree where our events live
  void listTrees(const char* options);

protected:

  // Uses the input event tree to discover the list of components
  Int_t buildComponents(TTree& eventTree);

  // Latch the values in a set of TTree
  Bool_t set(std::vector<TTree*>& trees);
 
  // Tree manipulation, called by openFile and makeFile
  Int_t makeBranches(TTree& eventTree, TTree& linkTree, TTree& fileTree, Int_t bufsize = 32000) const;
  Int_t attachToTree(TTree& eventTree, TTree& linkTree, TTree& fileTree);

  // cleanup in case there were problems with files
  void cleanup();

  // make sure that files exist.  pre-requisite for processing
  Bool_t checkDataFiles();

private:
  
  //disable copying and assignment
  CompositeEventList(const CompositeEventList& other);
  CompositeEventList& operator=(const CompositeEventList& other);

  // Data
  TTree*                                  _eventTree;  //! where the event data is
  TTree*                                  _linkTree;   //! where the links between event and meta data are
  TTree*                                  _fileTree;   //! where the list of files is
  
  DataHandle<Long64_t>                    _eventIndex; //! index of the current event in the _eventTree
  DataHandle<Long64_t>                    _metaIndex;  //! index of the current entry in the _fileTree
  DataHandle<Long64_t>                    _metaOffset; //! index of the first event with the current entry in _fileTree

  std::vector<TreeAndComponent>           _compList;   //! vector of components and associated TTrees
  std::vector<std::string >               _compNames;  //! names of components
  std::map<std::string,UInt_t>            _compMap;    //! components lookup map to get the index in the list above from the name

  ClassDef(CompositeEventList,0) // Base class for templates for handling simple types
};


#endif

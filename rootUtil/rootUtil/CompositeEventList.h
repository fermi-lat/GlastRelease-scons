
#ifndef CompositeEventList_h
#define CompositeEventList_h

/*
* Project: GLAST
* Package: rootUtil
* File: $Id$
* Authors:
*   EC, Eric Charles , SLAC, echarles@slac.stanford.edu
*   DC, David Chamont, LLR , chamont@poly.in2p3.fr
*
* Copyright (c) 2007
*                   Regents of Stanford University. All rights reserved.
*
*/

#include "CelEventInfo.h"
#include "DataHandle.h"
class CelEventComponent ;

#include <TObject.h>
#include <TString.h>
class TTree ;
class TFile ;
class TChain ;
class TVirtualIndex ;
class TCollection ;
class TObjArray ;

#include <map>
#include <vector>

//
// The class CompositeEventList permits to represent in memory
// a list of composite events, where each event is described as
// a set of entries in several separate trees. Each tree contains
// a different kind of data, and is also called a "component".
// Unlike using friend trees, this can describe sparse
// collections without requiring you to build a lookup index.
//
// An instance of CompositeEventList can be read from and/or
// written to a ROOT file. We should also consider the construction
// of a CompositeEventList on-the-fly, so to support the good old
// way to read/write events as was done by former RootIo writers
// and readers.
//
// The information is stored on 3 trees
// 
//   Events: 5 branches total, 1 entry per event
//      Event_Index/L     -> Index of the current composite event
//      Event_RunId/L     -> Unique? id of the corresponding run
//      Event_EventId/L   -> Unique? id of the corresponding event
//      Event_FileSetIndex/L  -> Index of the associated set of files and trees
//      Event_FileSetOffset/L -> Number of events in previous sets
// 
//   FileAndTreeSets: 5 Branches per input component, 1 entry per set
//      Comp_SetSize/s    -> Number of Files and Trees in the current set
//      Comp_FileNames    -> TClonesArray<TObjString> with the file names
//      Comp_TreeNames    -> TClonesArray<TObjString> with the tree names (should all be the same ?)
//      Comp_TreesSize/L  -> Total number of entries in all trees from the current set
//      Comp_TreeOffsets  -> TArrayL64 with the offsets for each Tree in the current set
//
//   ComponentEntries:  2 Branches per input component, 1 entry per event
//      Comp_EntryIndex/L  -> Index of the entry in the tree it lives on
//      Comp_TreeIndex/s -> Index of the tree in the associated set of files and trees
//
// Note : the couple (run_id/event_id) is expected to be unique for GLAST real data. For
// what concerns the simulation data, the triplet {task_name,run_id,event_id} should be
// unique.
//   

class CompositeEventList : public TObject
 {

  public :

    // constructors
    CompositeEventList() ; 
    CompositeEventList( TTree & eventTree, TTree & fileTree, TTree & entryTree ) ;
    virtual ~CompositeEventList() ;

    // Add a component by name.  This is only mandatory when writing.
    // On read these are discovered if not predeclared.
    UInt_t addComponent( const TString & name) ;
    
    // Make a new CEL ROOT file and its trees
    TFile * makeFile( const TString & fileName, const Char_t * options) ;
    
    // Open an existing CEL ROOT file
    TFile * openFile( const TString & fileName ) ;

    
    
    
  
  // Set up an event.  Grab the status of a set of TTrees
  Long64_t fillEvent(TObjArray& trees);
  Long64_t fillEvent(std::vector<TTree*>& trees);
  Long64_t fillFileAndTreeSet();
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
    // Number of components
    UInt_t numComponents() const { return _compList.size() ; }
    // Number of events in this cel
    Long64_t numEvents() const;
    // Get the index of the current event
    Long64_t eventIndex() const { return _currentEvent.eventIndex() ; }
    
    /// USEFUL IN PUBLIC ???
    // Get the index of the current entry in the file tree
    Long64_t fileSetIndex() const { return _currentEvent.fileSetIndex() ; }
    // Get the offset of the current entry in the file tree
    Long64_t fileSetOffset() const { return _currentEvent.fileSetOffset() ; }
    // Return the tree with the Event component entries
    TTree* entryTree() { return _entryTree; }
    // Return the tree with the links between events and meta data
    TTree* linkTree() { return _eventTree; }
    // Return the tree with the names of data files
    TTree* fileTree() { return _fileTree; }
    // Get the index of a component by name.  Returns 0xFFFFFFFF if no
    // such component exists
    UInt_t componentIndex( const TString & name )
     {
      std::map<TString,UInt_t>::const_iterator itrFind = _compMap.find(name) ;
      return itrFind != _compMap.end() ? itrFind->second : 0xFFFFFFFF ;
     } 
    // Get a Tree that is being read
    TTree * getTree(UInt_t index) const
     { return (index < _compList.size()) ? _compList[index].first : 0 ;  }
    // Get an event Component that is being read
    CelEventComponent * getComponent(UInt_t index) const
     { return (index < _compList.size()) ? _compList[index].second : 0 ; }

    
    
  // Printing
  // Dump a set of event component pointers and the list of TTree where they live
  void printout(const char* options, UInt_t nEvent=10, UInt_t startEvent=0);
  // Dump a single event
  void dumpEvent(const char* options);
  // Dump the list of TTree where our events live
  void listTrees(const char* options);

  
  protected :

    // Uses the input tree to discover the list of components
    Int_t buildComponents( TTree & entryTree ) ;

    // Latch the values in a set of TTree
    Bool_t set(std::vector<TTree*>& dataTrees) ;
 
    // Manipulation of cel internal trees
    void deleteCelTrees() ;
    Bool_t checkCelTrees() ;
    Int_t makeBranches( TTree & eventTree, TTree & fileTree, TTree & eventTree, Int_t bufsize = 32000) const;
    Int_t attachToTree( TTree & eventTree, TTree & fileTree, TTree & eventTree ) ;


  private :
  
	typedef std::pair<TTree*,CelEventComponent*> CelTreeAndComponent ;

    // disable copying and assignment
    CompositeEventList( const CompositeEventList & ) ;
    CompositeEventList & operator=( const CompositeEventList & ) ;

    // data
    TTree * _eventTree ; 
    TTree * _fileTree ; 
    TTree * _entryTree ;  
    CelEventInfo _currentEvent ;                      //! current link : event index, file-tree-name index
    std::vector<CelTreeAndComponent> _compList ;   //! vector of components and associated TTrees
    std::vector<TString>             _compNames ;  //! names of components
    std::map<TString,UInt_t>         _compMap ;    //! components lookup map to get the index in the list above from the name

    ClassDef(CompositeEventList,0)

 } ;


#endif

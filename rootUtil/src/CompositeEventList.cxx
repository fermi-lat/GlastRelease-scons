// -*- Mode: c++ -*-
#ifndef CompositeEventList_cxx
#define CompositeEventList_cxx
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
#include "rootUtil/CompositeEventList.h"

// c++/stl headers
#include <iostream>
#include <assert.h>

// ROOT Headers
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TObjArray.h>

// Other headers from this package
#include "rootUtil/TreeRefHandle.h"
#include "rootUtil/CelComponent.h"
#include "rootUtil/CelIndex.h"
#include "rootUtil/FileUtil.h"


ClassImp(CompositeEventList);


CompositeEventList::CompositeEventList():
  BranchGroup(),
  TObject(),
  _eventTree(0),
  _linkTree(0),
  _fileTree(0),
  _eventIndex(-1,*this,"EventIndex"),
  _metaIndex(0,*this,"MetaIndex"),
  _metaOffset(0,*this,"MetaOffset"){	
  // Default c'tor.  
}
	    
CompositeEventList::CompositeEventList(TTree& eventTree, TTree& linkTree, TTree& fileTree):
  BranchGroup(),
  TObject(),
  _eventTree(&eventTree),
  _linkTree(&linkTree),
  _fileTree(&fileTree),
  _eventIndex(-1,*this,"EventIndex"),
  _metaIndex(0,*this,"MetaIndex"),
  _metaOffset(0,*this,"MetaOffset"){	
  // Build a pointer skim from the three trees
  attachToTree(eventTree,linkTree,fileTree);
}

CompositeEventList::~CompositeEventList(){
  // Standard c'tor.  CelComponent Name is set
  cleanup();
}


TFile* CompositeEventList::makeFile(const Char_t* fileName, const Char_t* options) {
  // Make a new file.  Will also delare TTree for storing the pointers
  //
  // This makes two TTrees in the new file.
  //   Events -> One entry per event.  Has the indices to point to events.
  //   Files  -> One entry per input collection.  Stores the File names and TTree names. 
  //
  TDirectory* oldDir = gDirectory;
  TFile* newFile = TFile::Open(fileName,options);
  if ( 0 == newFile ) {
    // Root has already warned.  Just return NULL
    return 0;
  }
  newFile->cd();
  _eventTree = new TTree("Events","Events");
  _linkTree = new TTree("Links","Links");
  _fileTree = new TTree("Files","Files");
  oldDir->cd();
  Int_t check = makeBranches(*_eventTree,*_linkTree,*_fileTree,32000);
  if ( check < 0 ) {
    // Failed to make branches for some reason.  
    // Already warned.  clean up and return NULL
    cleanup();
    delete newFile;
    return 0;
  }
  return newFile;
}

TFile* CompositeEventList::openFile(const Char_t* fileName){
  // Open an existing file which contains a pointer skim
  //
  // This will warn and return NULL if 
  //   a) the file doesn't exist
  //   b) the file doesn't contain a pointer skim
  TDirectory* oldDir = gDirectory;
  TFile* newFile = TFile::Open(fileName);
  if ( 0 == newFile ) return 0;
  newFile->cd();
  _eventTree = dynamic_cast<TTree*>(newFile->Get("Events"));
  _linkTree = dynamic_cast<TTree*>(newFile->Get("Links"));
  _fileTree = dynamic_cast<TTree*>(newFile->Get("Files"));
  if  ( _eventTree == 0 || _linkTree == 0 || _fileTree == 0 ) {
    std::cerr << "Not TTrees" << std::endl;
    cleanup();
    delete newFile;
  }
  if ( _compList.size() == 0 ) {
    buildComponents(*_eventTree);
  }
  oldDir->cd();
  Int_t check = attachToTree(*_eventTree,*_linkTree,*_fileTree);
  if ( check < 0 ) {
    std::cerr << "File is not a valid skim" << std::endl;
    cleanup();
    delete newFile;   
    return 0;
  }
  return newFile;
}

UInt_t CompositeEventList::addComponent(const std::string& name){
  // Add a component by name.  This is only need when writing.  On read these are discovered.
  //
  // Returns the component index value
  CelComponent* comp = new CelComponent(name);
  TTree* t(0);
  TreeAndComponent tc(t,comp);
  // David: below, it should rather be _compList, no ?
  // (besides, it works this way, lucky guys)
  UInt_t retValue = _compMap.size();
  _compMap[name] = retValue;
  _compNames.push_back(name);
  _compList.push_back(tc);
  return retValue;
}

Long64_t CompositeEventList::fillEvent(TObjArray& trees) {
  // Set up an event.  Grab the status of a set of TTrees
  //
  // Returns the entry number of the event that has just been written
  // Returns -1 if it failed to fill the event
  // This is here since CINT screws up with std::vector<TTree*> 
  std::vector<TTree*> v;
  for ( Int_t i(0); i < trees.GetEntries(); i++ ) {
    TTree* t = static_cast<TTree*>(trees.UncheckedAt(i));
    v.push_back(t);
  }
  return fillEvent(v);
}

Long64_t CompositeEventList::fillEvent(std::vector<TTree*>& trees) {
  // Set up an event.  Grab the status of a set of TTrees
  //
  // Returns the entry number of the event that has just been written  
  if ( ! checkDataFiles() ) return -1;
  if ( ! set(trees) ) return -1;
  _eventIndex = _eventIndex + 1;
  _eventTree->Fill();
  _linkTree->Fill();
  return _eventIndex;
}

Long64_t CompositeEventList::fillMeta() {
  // Store up the List of TTrees that have been used so far
  //
  // Returns the entry number of the entry that has just been written    
  if ( ! checkDataFiles() ) return -1;
  _fileTree->Fill();
  Long64_t retValue = _metaIndex;
  _metaIndex = _metaIndex + 1;
  _metaOffset = _eventIndex;
  return retValue;
}


Int_t CompositeEventList::read(Long64_t iEvt) {
  // Read an event by Index
  // 
  // If successful, returns the number of bytes read
  // Failure codes:
  //   -1 -> Input Trees do not exist.  Probably not yet attached
  //   -2 -> Can't read entry on Event input tree
  //   -3 -> Can't read entry on Link input tree
  //   -4 -> Can't read entry on File input tree
  //   -5 -> Can't find a component
  //   -6 -> Can't read a component input tree
  Int_t total = readEventAddress(iEvt);
  if ( total < 0 ) return total;

  for ( std::vector< TreeAndComponent >::iterator itr = _compList.begin();
	itr != _compList.end(); itr++ ) {
    CelComponent* comp = itr->second;
    if ( comp == 0 ) return -5;
    Int_t nRead = comp->read();
    if ( nRead < 0 ) return -6;
    TTree* t = comp->getTree();
    itr->first = t;
    total += nRead;
  }
  return total;
}


Int_t CompositeEventList::readEventAddress(Long64_t iEvt) {
  // read only the event tree data
  //
  // If successful, returns the number of bytes read
  // Failure codes:
  //   -1 -> Input Trees do not exist.  Probably not yet attached
  //   -2 -> Can't read entry on Event input tree
  //   -3 -> Can't read entry on Link input tree
  //   -4 -> Can't read entry on File input tree  
  if ( !checkDataFiles() ) return -1;

  if ( iEvt == _eventIndex ) {
    // event already read,
    return 0;
  }

  if ( iEvt == 0 ) {
    ;
  } else if ( iEvt % 10000 == 0 ) {
    std::cout << 'x' << std::flush;
  } else if ( iEvt % 1000 == 0 ) {
    std::cout << '.' << std::flush;
  } 
  _eventIndex = iEvt;
  Int_t total = _eventTree->GetEntry(iEvt);
  if ( total < 0 ) return -2;
  Int_t link = _linkTree->GetEntry(iEvt);
  if ( link < 0 ) return -3;
  total += link;

  Long64_t mEvt = _metaIndex;

  Int_t nMeta = _fileTree->GetEntry(mEvt);
  if ( nMeta < 0 ) return -4;
  total += nMeta;

  return total;
}


TChain* CompositeEventList::buildChain(UInt_t index) {
  // Build the TChain for a component. 
  if ( ! checkDataFiles() ) {
    std::cerr << "No Data Files" << std::endl;
    return 0;
  }
  CelComponent* comp = getComponent(index);
  if ( 0 == comp ) {
    std::cerr << "No component with index " << index << std::endl;
    return 0;
  }
  TChain* chain(0);
  std::string chainName;
  for ( Long64_t i(0); i < _fileTree->GetEntries(); i++ ) {
    Int_t nB = _fileTree->GetEntry(i);    
    if ( nB < 0 ) {
      std::cerr << "Could not read Tree entry " << i << " for component " << _compNames[i] << std::endl;
      delete chain;
      return 0;
    }
    // this will build a new chain if needed
    Bool_t isOk = comp->addToChain(chain);
    if ( !isOk ) {
      std::cerr << "Could not add tree " << i << " to chain for component " << _compNames[i] << std::endl;
      delete chain;
      return 0;
    }
  }
  return chain;
}

TVirtualIndex* CompositeEventList::buildEventIndex(UInt_t index, Long64_t& offset, TTree* tree) {
  // Build the TChain for a component. 
  if ( ! checkDataFiles() ) {
    // FIXME: Warn?
    return 0;
  }
  TVirtualIndex* vIdx = CelIndex::buildIndex(*this,_compNames[index],tree,offset);
  return vIdx;
}

TChain* CompositeEventList::buildLinks(TObjArray* chainList, Bool_t setFriends){
  // Build the event chain

  if ( 0 == _eventTree ) return 0;
  TFile* f = FileUtil::getFile(*_eventTree);
  if ( 0 == f ) return 0;
    
  TChain* newChain = new TChain("Events");  
  Int_t check = newChain->Add(f->GetName());
  if ( check < 0 ) {
    delete newChain;
    return 0;
  }

  // Build all the Chains
  Long64_t offset(0);
  for ( UInt_t i(0); i < nComponent(); i++ ) {
    std::cout << "Building links for " << _compNames[i] << std::endl;
    TChain* c = buildChain(i);
    c->SetBranchStatus("*",0);
    if ( 0 == c ) {
      // FIXME: warn?
      return 0;
    }
    TVirtualIndex* vIdx = buildEventIndex(i,offset,c);
    if ( 0 == vIdx ) {
      // FIXME: warn?
      return 0;
    }     
    if ( 0 != chainList ) {
      chainList->AddLast(c);
    }
    if ( setFriends ) {
       newChain->AddFriend(c,_compNames[i].c_str());
    }
  }
  return newChain;
}


Long64_t CompositeEventList::entries() const { 
  // return the number of entries
  return _eventTree != 0 ? _eventTree->GetEntries() : 0;
}

void CompositeEventList::printout(const char* options, UInt_t nEvent, UInt_t startEvent){
  // Dump a set of event component pointers and the list of TTree where they live
  //
  // One line per event
  UInt_t lastEvt = startEvent+nEvent;
  std::cout << std::endl << "Printing Events: " << std::endl;
  for ( UInt_t iEvt(startEvent); iEvt < lastEvt; iEvt++ ) {
    read(iEvt);
    dumpEvent(options);
  }

  std::cout << std::endl << "Printing Trees: " << std::endl;
  listTrees(options);
}

void CompositeEventList::dumpEvent(const char* options){
  // Dump a single event on one line
  // 
  // The format is:
  //  Event:     comp1 tree0:event0   comp1 tree1:event1 ... 
  std::cout << "Event " << _eventIndex << ":\t";
  for ( std::vector< TreeAndComponent >::const_iterator itr = _compList.begin();
	itr != _compList.end(); itr++ ) {
    const CelComponent* comp = itr->second;
    if ( comp == 0 ) return;
    std::cout << comp->componentName() << ' ';
    comp->dumpEvent(options);
    std::cout << "\t";
  } 
  std::cout << std::endl;
}

void CompositeEventList::listTrees(const char* options){
  // Dump the list of TTree where our events live
  //
  // One tree per line, the format is:
  //   File Tree
  for ( std::vector< TreeAndComponent >::const_iterator itr = _compList.begin();
	itr != _compList.end(); itr++ ) {
    const CelComponent* comp = itr->second;
    if ( comp == 0 ) return;
    std::cout << comp->componentName() << std::endl;
    comp->listTrees(options);
  }   
}



Int_t CompositeEventList::buildComponents(TTree& tree){
  // Uses the input event tree to discover the list of components
  //
  // Called by openFile()
  TObjArray* array = tree.GetListOfBranches();
  for ( Int_t i(0); i < array->GetEntries(); i++ ) {
    TObject* obj = array->UncheckedAt(i);
    std::string name(obj->GetName());    
    std::string::size_type find =  name.find("EvtIndex");
    if ( find == name.npos ) continue;
    std::string compName(name,0,find);
    addComponent(compName);
  }
  return nComponent();
}


Bool_t CompositeEventList::set(std::vector<TTree*>& trees) {
  // Latch the values in a set of TTree
  //
  // Called by fillEvent()
  assert ( checkDataFiles() ) ;
  if (trees.size() != _compList.size()) {
    std::cerr << "Wrong number of trees.  " 
	      << _compList.size() << " components and " << trees.size() << " trees." << std::endl;
    return kFALSE;
  }
  
  UInt_t idx(0);  
  for ( std::vector< TTree* >::iterator itr = trees.begin();
	itr != trees.end(); itr++ ) {  
    CelComponent* comp = getComponent(idx);
    assert ( 0 != comp );
    comp->set(**itr);
    idx++;
  }
  return kTRUE;
}


Int_t CompositeEventList::makeBranches(TTree& eventTree, TTree& linkTree,  TTree& fileTree, Int_t bufsize) const {
  // Make the Branches we store everything on
  // 
  // Called by makeFile()
  Int_t total = BranchGroup::makeBranches(linkTree,"Link_",bufsize);
  if ( total < 0 ) return total;
  for ( std::vector< TreeAndComponent >::const_iterator itr = _compList.begin();
	itr != _compList.end(); itr++ ) {
    const CelComponent* comp = itr->second;
    if ( comp == 0 ) return -1;
    Int_t nMake = comp->makeBranches(eventTree,fileTree,bufsize);
    if ( nMake < 0 ) return nMake;
    total += nMake;
  }
  return total;
}



Int_t CompositeEventList::attachToTree(TTree& eventTree, TTree& linkTree,  TTree& fileTree) {
  // Attach to the branches where all the infomation is stored
  // 
  // Called by openFile()
  Int_t total = BranchGroup::attachToTree(linkTree,"Link_");
  if ( total < 0 ) {
    std::cerr << "Failed to open Link Tree" << std::endl;    
    return total;
  }
  UInt_t idx(0);
  for ( std::vector< TreeAndComponent >::const_iterator itr = _compList.begin();
	itr != _compList.end(); itr++ ) {
    CelComponent* comp = itr->second;
    if ( comp == 0 ) { 
      std::cerr << "Failed to get component " << idx << std::endl;
      return -1;
    }
    Int_t nB = comp->attachToTree(eventTree,fileTree);
    if ( nB < 0 ) {
      std::cerr << "Failed to attached component " << comp->componentName() << " to trees " 
		<< eventTree.GetName() << ' ' << fileTree.GetName() << std::endl;
      return nB;
    }
    total += nB;
    idx++;
  }
  return total;
}

void CompositeEventList::cleanup() {
  // cleanup in case there were problems with files
  delete _eventTree; _eventTree = 0;
  delete _linkTree; _linkTree = 0;
  delete _fileTree; _fileTree = 0;
}

Bool_t CompositeEventList::checkDataFiles() {
  // make sure that files exist.  pre-requisite for processing
  if ( 0 == _eventTree ) return kFALSE;
  if ( 0 == _linkTree ) return kFALSE;  
  if ( 0 == _fileTree ) return kFALSE;  
  return kTRUE;
}

#endif

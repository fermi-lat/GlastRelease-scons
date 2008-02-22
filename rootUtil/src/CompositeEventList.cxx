
/*
* Project: GLAST
* Package: rootUtil
* File: CompositeEventList.cxx
* Authors:
*   DC, David Chamont, LLR,  chamont@llr.in2p3.fr
*   EC, Eric Charles,  SLAC, echarles@slac.stanford.edu
*
*/

#include "rootUtil/CompositeEventList.h"
#include "rootUtil/CelEventComponent.h"
#include "rootUtil/CelIndex.h"
#include "rootUtil/RuChain.h"
#include "rootUtil/RuUtil.h"

#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TSystem.h>
#include <Riostream.h>

#include <assert.h>



//====================================================================
// globals and low level utilities
//====================================================================


ClassImp(CompositeEventList) ;

static std::string __entryTreeName("EventEntries") ;
static std::string __linkTreeName("EventLinks") ;
static std::string __fileTreeName("FileAndTreeSets") ;
static std::string __offsetTreeName("FileAndTreeOffsets") ;

const UInt_t CompositeEventList::COMPONENT_UNDEFINED = static_cast<UInt_t>(-1) ;

Bool_t CompositeEventList::checkCelTree
 ( const TString & caller, TTree * tree, const TString & name, Bool_t error )
 {
  if (tree==0)
   {
    if (error==kTRUE)
     {
      std::cerr<<"["<<caller<<"] "<<"lacking tree "<<name<<std::endl ;
      _isOk = kFALSE ;
      deleteCurrentFile() ;
     }
    else
     {
      std::cout<<"["<<caller<<"] "<<"no tree "<<name<<std::endl ;
     }
    return kFALSE ;
   }
  if (name!=tree->GetName())
   {
    std::cerr<<"["<<caller<<"] "
      << "a tree is called "<<tree->GetName()
      << " instead of "<<name
      << std::endl ;
    _isOk = kFALSE ;
    deleteCurrentFile() ;
    return kFALSE ;
   }
  return kTRUE ;
 }

Bool_t CompositeEventList::checkCelOk( const TString & caller ) const
 {
  if (_isOk==kFALSE)
   { std::cerr<<"["<<caller<<"] cel not OK"<<std::endl ; }
  return _isOk ;
 }

Bool_t CompositeEventList::checkCelPrepared( const TString & caller ) const
 {
  if (checkCelOk(caller)==kFALSE)
   { return kFALSE ; }
  
  if (_isPrepared==kFALSE)
   {
    if (_openingOptions=="RECREATE")
     {
      std::cerr
        << "[CompositeEventList::checkCelPrepared] cannot prepare a const CEL for "
        << _openingOptions
        << std::endl ;
      return kFALSE ;
     } 
    else if (_openingOptions=="READ")
     { prepareRead() ; }
    else
     {
      std::cerr
        << "[CompositeEventList::checkCelPrepared] no preparation known for "
        << _openingOptions
        << std::endl ;
      return kFALSE ;
     }
   }
  if (_isPrepared==kFALSE)
   {
    _isOk = kFALSE ;
    std::cerr<<"["<<caller<<"] cel not prepared"<<std::endl ;
    return kFALSE ;
   }
  return _isPrepared ;
 }

Bool_t CompositeEventList::checkCelPrepared( const TString & caller )
 {
  if (checkCelOk(caller)==kFALSE)
   { return kFALSE ; }
  
  if (_isPrepared==kTRUE)
   { return kTRUE ; }
  if (_openingOptions!="RECREATE")
   { return (static_cast<const CompositeEventList *>(this))->checkCelPrepared(caller) ; }
  prepareRecreate() ;
  if (_isPrepared==kFALSE)
   {
    _isOk = kFALSE ;
    std::cerr<<"["<<caller<<"] CEL not prepared"<<std::endl ;
    return kFALSE ;
   }
  return _isPrepared ;
 }

Bool_t CompositeEventList::checkCelTrees( const TString & caller )
 {
  if (_currentFile==0)
   {
    _isOk = kFALSE ;
    return kFALSE ;
   }
  if (checkCelTree(caller,_entryTree,__entryTreeName,kTRUE)==kFALSE) return kFALSE ;
  if (checkCelTree(caller,_linkTree,__linkTreeName,kTRUE)==kFALSE) return kFALSE ;
  if (checkCelTree(caller,_fileTree,__fileTreeName,kTRUE)==kFALSE) return kFALSE ;
  checkCelTree(caller,_offsetTree,__offsetTreeName,kFALSE) ;
  return kTRUE ;
 }


 
//====================================================================
// construction
//====================================================================


CompositeEventList::CompositeEventList
 ( const TString & celFileName,
   const TString & options,
   const TObjArray * componentNames )
 : _isOk(kTRUE), _isPrepared(kFALSE), _fileName(celFileName), _openingOptions(options),
   _currentFile(0), _entryTree(0), _linkTree(0), _fileTree(0), _offsetTree(0)
 {
  TString currentMethod("CompositeEventList::CompositeEventList") ;
  
  // transient not implemented
  if (celFileName=="")
   {
    std::cerr<<"["<<currentMethod<<"] "
      <<"transient CELs not yet implemented "
      <<std::endl ;
    _isOk = kFALSE ;
   }
  
  // check option
  if ( (_openingOptions != "RECREATE") &&
	   (_openingOptions != "READ") )
   {
    _openingOptions = "" ;
    std::cerr<<"["<<currentMethod<<"] "
      << "the option " << options << " is not yet implemented"
      << std::endl ;
    _isOk = kFALSE ;
   }
  
  // expand file name
  TString fileName = celFileName ;
  if (gSystem->ExpandPathName(fileName)==kTRUE)
   {
    std::cerr<<"["<<currentMethod<<"] "
      << "Failed to expand some variable in filename "
      << celFileName
      << std::endl ;
    _isOk = kFALSE ;
   }
  
  // open file
  rootUtil::AutoCd cdBackup ; // will restore old global directory when deleted
  _currentFile = TFile::Open(fileName,options) ;
  if ((_currentFile==0)||(!_currentFile->IsOpen()))
   {
    std::cerr<<"["<<currentMethod<<"] "
      << "Failed to open file "<< fileName
      << " in "<< _openingOptions << " mode"
      << std::endl ;
	delete _currentFile ;
	_currentFile = 0 ;
	_isOk = kFALSE ; ;
   } 
  
  // prepare all that can be before all the components
  // are known : the trees.
  if (_openingOptions=="RECREATE")
   {
    // make trees
    _entryTree = new TTree(__entryTreeName.c_str(),__entryTreeName.c_str()) ;
    _linkTree = new TTree(__linkTreeName.c_str(),__linkTreeName.c_str()) ;
    _fileTree = new TTree(__fileTreeName.c_str(),__fileTreeName.c_str()) ;
    _offsetTree = new TTree(__offsetTreeName.c_str(),__offsetTreeName.c_str()) ;
   }
  else
   {
    // explore
    _entryTree = dynamic_cast<TTree*>(_currentFile->Get(__entryTreeName.c_str())) ;
    _linkTree = dynamic_cast<TTree*>(_currentFile->Get(__linkTreeName.c_str())) ;
    _fileTree = dynamic_cast<TTree*>(_currentFile->Get(__fileTreeName.c_str())) ;
    _offsetTree = dynamic_cast<TTree*>(_currentFile->Get(__offsetTreeName.c_str())) ;
   }
  checkCelTrees(currentMethod) ;
  
  // some or all the components names can be given
  // to the constructors
  if (componentNames!=0)
   { declareComponents(componentNames) ; }

 }
	    
void CompositeEventList::deleteCurrentFile()
 {
  delete _currentFile ;
  _currentFile = 0 ;
  _entryTree = 0  ;
  _linkTree = 0   ;
  _fileTree = 0   ;
  _offsetTree = 0 ;
  _openingOptions = "" ;
 }

Bool_t CompositeEventList::isOk()
 { return _isOk ; }

CompositeEventList::~CompositeEventList()
 {
  // no file open
  if (_currentFile==0)
   { return ; }
  
  // eventually needs writing
  if (_openingOptions!="READ")
   { writeAndClose() ; }
  else
   { deleteCurrentFile() ; }
 }



//====================================================================
// Accessors
//====================================================================


UInt_t CompositeEventList::numComponents() const
 { return _compList.size() ; }

const TString & CompositeEventList::componentName( UInt_t componentIndex ) const
 { return _compNames[componentIndex] ; }

UInt_t CompositeEventList::componentIndex( const TString & componentName ) const
 {
  std::map<TString,UInt_t>::const_iterator itrFind = _compMap.find(componentName) ;
  return itrFind != _compMap.end() ? itrFind->second : COMPONENT_UNDEFINED ;
 }

Long64_t CompositeEventList::numEvents() const
 { return (_linkTree!=0?_linkTree->GetEntries():0) ; }

Long64_t CompositeEventList::numFileAndTreeSets() const
 { return (_fileTree!=0?_fileTree->GetEntries():0) ; }

const TString & CompositeEventList::fileName() const
 { return _fileName ; }

Long64_t CompositeEventList::currentEventIndex() const
 { return _currentLink.eventIndex() ; }

Long64_t CompositeEventList::currentSetIndex() const
 { return _currentLink.setIndex() ; }



//====================================================================
// components setup
//====================================================================


UInt_t CompositeEventList::declareComponent( const TString & name )
 {
  // Add a component by name.  This is only need when writing.
  // On read these are discovered.
  // Returns the component index value*
//  std::cout<<"[CompositeEventList::declareComponent] declaring component "<<name<<std::endl ;
  UInt_t retValue = _compList.size() ;
  _compList.push_back(new CelEventComponent(name)) ;
  _compNames.push_back(name) ;
  _compMap[name] = retValue ;
  return retValue ;
 }

UInt_t CompositeEventList::declareComponents( const TObjArray * names )
 {
  Int_t i ;
  for ( i=0 ; i < names->GetEntries(); i++ )
   {
    TObjString * s = static_cast<TObjString*>(names->UncheckedAt(i)) ;
    declareComponent(s->GetString()) ;
   }
  return numComponents() ;
 }

// Uses the input file tree to discover the list of components
UInt_t CompositeEventList::discoverComponents()
 {
  TObjArray * array = _fileTree->GetListOfBranches() ;
  Int_t i ;
  for ( i = 0 ; i < array->GetEntries() ; i++ )
   {
    TObject * obj = array->UncheckedAt(i) ;
    TString name = obj->GetName() ;
    Ssiz_t find = name.Index("_Set_Size") ;
    if (find==kNPOS) continue ;
    name.Remove(find) ;
    declareComponent(name) ; 
    //std::string name(obj->GetName()) ;    
    //std::string::size_type find =  name.find("EntryIndex") ;    
    //if ( find == name.npos ) continue ;   
    //std::string compName(name,0,find) ; 
    //addComponent(compName) ;
   }
  return numComponents() ;
 }

CelEventComponent * CompositeEventList::getComponent( UInt_t index ) const
 { return (index < _compList.size()) ? _compList[index] : 0 ; }

CelEventComponent * CompositeEventList::getComponent( const TString & name ) const
 {
  UInt_t index = componentIndex(name) ;
  return index != COMPONENT_UNDEFINED ? getComponent(index) : 0 ;
 }

TTree * CompositeEventList::getTree( UInt_t index ) const
 { return (index < _compList.size()) ? _compList[index]->getTree() : 0 ; }

TTree * CompositeEventList::getTree( const TString & name ) const
 {
  UInt_t index = componentIndex(name) ;
  return index != COMPONENT_UNDEFINED ? getTree(index) : 0 ;
 }

//void CompositeEventList::prepare()
// {
//  std::cout<<"[CompositeEventList::prepare] begin"<<std::endl ;
//
//  // check components
//   
//  // delegate real preparation
//  if (_openingOptions=="RECREATE")
//   { prepareRecreate() ; } 
//  else if (_openingOptions=="READ")
//   { prepareRead() ; }
//  else
//   {
//    std::cerr
//      << "[CompositeEventList::prepare] no preparation known for "
//      << _openingOptions
//      << std::endl ;
//   }
// }


//====================================================================
// write use-case
//====================================================================


void CompositeEventList::prepareRecreate()
 {
  // check status
  TString currentMethod("CompositeEventList::prepareRecreate") ;
  if (_isPrepared==kTRUE)
   {
    std::cout
      << "["<<currentMethod<<"] "
      << _fileName << "already prepared ?! "<< std::endl ;
    return ;
   }
        
  // check opening mode
  if (_openingOptions!="RECREATE")
   {
    std::cerr
      << "[CompositeEventList::prepareRecreate] "
      << "inconsistent with openingOptions "<< _openingOptions
      << std::endl ;
    deleteCurrentFile() ;
    _isOk = kFALSE ;
   }
        
  // make branches
  Int_t check = makeCelBranches(_entryTree,_linkTree,_fileTree,_offsetTree,32000) ;
  if ( check < 0 )
   {
    // Failed to make branches for some reason.  
    // Already warned.  clean up and return NULL
    deleteCurrentFile() ;
    _isOk = kFALSE ;
   }
  
  //
  if (_isOk==kFALSE)
   {
    std::cerr
      << "[CompositeEventList::prepareRecreate] preparation failed"
      << std::endl ;
   } 
  else
   {
    _isPrepared = kTRUE ;
   }
 }

Int_t CompositeEventList::makeCelBranches
 ( TTree * entryTree, TTree * linkTree,  TTree * fileTree,  TTree * offsetTree, Int_t bufsize ) const
 {
  Int_t n, total = 0 ;
  if (linkTree!=0)
   {
	n = _currentLink.makeBranches(*linkTree,0,bufsize) ;
	if ( n < 0 ) return n ;
    total += n ;
   }
  std::vector<CelEventComponent *>::const_iterator compItr ;
  for ( compItr = _compList.begin() ; compItr != _compList.end() ; compItr++ )
   {
    const CelEventComponent * comp = *compItr ;
    if ( comp == 0 ) return -1 ;
    n = comp->makeBranches(entryTree,fileTree,offsetTree,bufsize) ;
    if ( n < 0 ) return n ;
    total += n ;
   }
  return total ;
 }

Long64_t CompositeEventList::fillEvent( const TObjArray & trees )
 {
  std::vector<TTree*> v;
  for ( Int_t i(0); i < trees.GetEntries(); i++ )
   {
    TTree* t = static_cast<TTree*>(trees.UncheckedAt(i));
    v.push_back(t);
   }
  return fillEvent(v);
 }

Long64_t CompositeEventList::fillEvent( const std::vector<TChain *> & chains )
 {
  std::vector<TTree *> v ;
  std::vector<TTree *>::size_type i ;
  for ( i = 0 ; i < chains.size() ; i++ )
   { v.push_back(chains[i]->GetTree()) ; }
  return fillEvent(v) ;
 }

// Insert an event in the CEL. Grab the status of a set of TTrees
// Returns the entry number of the event that has just been written  
// Returns -1 if it failed to fill the event
// This is here since CINT screws up with std::vector<TTree*> 
Long64_t CompositeEventList::fillEvent( const std::vector<TTree *> & trees )
 {
  TString caller ("[CompositeEventList::fillEvent] ") ;

  if (trees.size() != _compList.size())
   {
    std::cerr<<caller
      << "Do not know how to match " << _compList.size() << " components "
      << "with " << trees.size() << " trees."
      << std::endl ;
    return -1 ;
   }
  
  // TODO : check tree and component names ?
  
  UInt_t idx(0) ;
  std::vector< TTree* >::const_iterator itr ;
  for ( itr = trees.begin(), idx = 0 ; itr != trees.end() ; itr++, idx++ )
   { if (fillEntry(idx,*itr)==kFALSE) return -1 ; }
  
  return fillEvent() ;
 }

Bool_t CompositeEventList::fillEntry( const TString & componentName, TTree * tree )
 { return fillEntry(componentIndex(componentName),tree) ; }

  // Register the current entry for a given component
// After one has called fillEntry() for all components,
// he must call fillEvent().
Bool_t CompositeEventList::fillEntry(  UInt_t componentIndex, TTree * tree )
 {
  TString caller ("CompositeEventList::fillEntry") ;
  if ( ! checkCelPrepared(caller) ) return kFALSE ;
    
  // TODO : check tree and component names ?
  CelEventComponent * comp = getComponent(componentIndex) ;
  assert ( 0 != comp ) ;
  comp->registerEntry(*tree) ;
  
  return kTRUE ;
}

// Insert an event in the CEL, after some calls to fillEntry().
// Returns the entry number of the event that has just been written  
// Returns -1 if it failed to fill the event
// This is here since CINT screws up with std::vector<TTree*> 
Long64_t CompositeEventList::fillEvent()
 {
  TString caller ("CompositeEventList::fillEvent") ;
  if ( ! checkCelPrepared(caller) ) return -1 ;
  
  //std::cout<<"["<<caller<<"] fill event "
  //  <<(_currentLink.eventIndex()+1)
  //  <<" with "<<trees.size()<<" components"
  //  <<std::endl ;

  if (_currentFile->TestBits(TFile::kWriteError))
   {
    std::cerr
      << "[CompositeEventList::fillEvent] kWriteError"
      << std::endl ;
    return -1 ;
   }
        
  _currentLink.incrementEventIndex() ;
  _entryTree->Fill() ;
  _linkTree->Fill() ;
  return _currentLink.eventIndex() ;
}

// Store up the List of TTrees that have been used so far
// Returns the entry number of the entry that has just been written    
Long64_t CompositeEventList::fillFileAndTreeSet()
 {
  TString caller ("CompositeEventList::fillFileAndTreeSet") ;
  if ( ! checkCelPrepared(caller) ) return -1 ;
  
  //std::cout<<"["<<caller<<"] begin"<<std::endl ;

  _fileTree->Fill() ;
  _offsetTree->Fill() ;

  Long64_t retValue = _currentLink.setIndex() ;
  _currentLink.incrementFileSetIndex() ;
  
  std::vector<CelEventComponent*>::iterator itr ;
  for ( itr = _compList.begin() ; itr != _compList.end() ; itr++ )
   {
    CelEventComponent * comp = *itr ;
    if ( comp == 0 ) return -1 ;
    comp->nextSet() ;
   }
  
  return retValue ;
}

void CompositeEventList::writeAndClose()
 {
  if (_openingOptions=="READ")
   {
	std::cerr
	  << "[CompositeEventList::writeAndClose] "
      << "Why ask to write a read file ??"
      << std::endl ;
   }
  
  if (_currentFile==0)
   {
	std::cerr
	  << "[CompositeEventList::writeAndClose] "
      << "Probaly already closed"
      << std::endl ;
   }
  
  TString caller("CompositeEventList::writeAndClose") ;
  if ( ! checkCelOk(caller) ) return ;
  if ( ! checkCelPrepared(caller) ) return ;
  
  std::cout<<"["<<caller<<"] begin"<<std::endl ;

  _currentFile->Write() ;
  deleteCurrentFile() ;
  _isOk = kTRUE ;
 }


 
//====================================================================
// read use-case
//====================================================================


void CompositeEventList::prepareRead() const
 {
  // check status
  if (_isPrepared==kTRUE)
   {
    std::cout
      << "[CompositeEventList::prepareRead] "
      << _fileName << "already prepared ?! "<< std::endl ;
    return ;
   }
        
  // check opening mode
  if (_openingOptions!="READ")
   {
    std::cerr
      << "[CompositeEventList::prepareRead] "
      << "inconsistent with openingOptions "<< _openingOptions
      << std::endl ;
    _isOk = kFALSE ;
    return ;
   }
  
  // components
  if (_compList.size()==0)
   { const_cast<CompositeEventList *>(this)->discoverComponents() ; }
	
  // attach
  Int_t check = attachToTree(_entryTree,_linkTree,_fileTree,_offsetTree) ;
  if (check<0)
   {
    std::cerr
	  << "[CompositeEventList::prepareRead] "
      << "File is not a valid composite event list"
      << std::endl ;
    _isOk = kFALSE ;
   }
  
  //
  if (_isOk==kFALSE)
   {
    std::cerr
      << "[CompositeEventList::prepareRead] preparation failed"
      << std::endl ;
   } 
  else
   {
    _isPrepared = kTRUE ;
   }
 }

Int_t CompositeEventList::attachToTree( TTree * entryTree, TTree * linkTree,  TTree * fileTree,  TTree * offsetTree ) const
 {
  // Attach to the branches where all the infomation is stored
  // 
  // Called by openFile()
  Int_t n, total = 0 ;
  if (linkTree!=0)
   {
    n = _currentLink.attachToTree(*linkTree) ;
    if ( n < 0 )
     {
      std::cerr << "Failed to attach Link Tree" << std::endl ;    
      return n ;
     }
    total += n ;
   }
  UInt_t idx = 0 ;
  std::vector<CelEventComponent*>::const_iterator itr ;
  for ( itr = _compList.begin() ; itr != _compList.end(); itr++ )
   {
    CelEventComponent * comp = *itr ;
    if ( comp == 0 )
     { 
      std::cerr << "Failed to get component " << idx << std::endl ;
      return -1 ;
     }
    Int_t nB = comp->attachToTree(entryTree,fileTree,offsetTree) ;
    if ( nB < 0 )
     {
      std::cerr
        << "Failed to attached component " << comp->componentName() << " to trees " 
		<< entryTree->GetName() << ' ' 
		<< fileTree->GetName() << ' ' 
		<< offsetTree->GetName() << std::endl ;
      return nB;
     }
    total += nB ;
    idx++;
   }
  return total ;
 }

// Build the TChain for a component. 
RuChain * CompositeEventList::newChain( UInt_t componentIndex ) const
 {
  TString caller ("CompositeEventList::newChain") ;
  if ( ! checkCelPrepared(caller) ) return 0 ;
  
  // get component
  CelEventComponent * comp = getComponent(componentIndex) ;
  if (comp==0)
   {
    std::cerr<<"No component with index "<<componentIndex<<std::endl ;
    return 0 ;
   }
  
  // the chain is made on the fly because we
  // do not know the tree name.
  RuChain * chain = 0 ;
  Long64_t i ;
  for ( i=0 ; i<_fileTree->GetEntries() ; i++ )
   {
    Int_t nB = _fileTree->GetEntry(i) ;    
    if (nB<0)
     {
      std::cerr
        << "Could not read _fileTree set " << i
        << " for component " << _compNames[i]
        << std::endl ;
      delete chain ;
      return 0 ;
     }
    // this will build a new chain if needed
    Bool_t isOk = comp->addSetToChain(chain) ;
    if (!isOk)
     {
      std::cerr
        << "Could not add tree " << i 
        << " to chain for component " << _compNames[i] 
        << std::endl ;
      delete chain ;
      return 0 ;
     }
   }
  return chain ;
 }

// If successful, returns the number of bytes read
// Failure codes:
//   -1 -> Al least one cel internal tree does not exist.
//   -2 -> Can't read entry on ComponentEntries
//   -3 -> Can't read entry on CompositeEvents
//   -4 -> Can't read entry on FileAndTreeSets
Int_t CompositeEventList::shallowRead( Long64_t eventIndex )
 {
  TString caller ("CompositeEventList::shallowRead") ;
  if ( ! checkCelOk(caller) ) return -1 ;
  if ( ! checkCelPrepared(caller) ) return -1 ;
    
  if ( eventIndex == _currentLink.eventIndex() )
   { return 0 ; }

  // is it useful below, for a shallow (=>fast) read ???
  if ( eventIndex != 0 )
   {
    if ( eventIndex % 10000 == 0 )
     { std::cout << 'x' << std::flush ; }
    else if ( eventIndex % 1000 == 0 )
     { std::cout << '.' << std::flush ; } 
   }
  
  //_currentLink.changeEventIndex(eventIndex) ; // ??
  
  // event entries
  Int_t totalRead = _entryTree->GetEntry(eventIndex) ;
  if ( totalRead < 0 ) return -2 ;
  Int_t eventInfoRead = _linkTree->GetEntry(eventIndex) ;
  if ( eventInfoRead < 0 ) return -3 ;
  totalRead += eventInfoRead ;

  Long64_t setIndex = _currentLink.setIndex() ;
  if (_fileTree->GetReadEntry()!=setIndex)
   {
    Int_t setRead = _fileTree->GetEntry(setIndex) ;
    if ( setRead < 0 ) return -4 ;
    totalRead += setRead ;
    Int_t offsetRead = _offsetTree->GetEntry(setIndex) ;
    if ( offsetRead < 0 ) return -4 ;
    totalRead += offsetRead ;
   }

  return totalRead ;
 }

Long64_t CompositeEventList::entryIndex( UInt_t componentIndex ) const
 {
  TString caller ("[CompositeEventList::entryIndex]") ;
  if ( ! checkCelOk(caller) ) return -1 ;
  if ( ! checkCelPrepared(caller) ) return -1 ;
  if (componentIndex==COMPONENT_UNDEFINED)
   {
    std::cerr<<caller<<" undefined component"<<std::endl ;
    return -1 ;
   }  
  CelEventComponent * comp = getComponent(componentIndex) ;
  return comp->currentIndexInChain() ;
 }

const TObjString * CompositeEventList::fileName( UInt_t componentIndex ) const
 {
  TString caller ("[CompositeEventList::fileName]") ;
  if ( ! checkCelOk(caller) ) return 0 ;
  if ( ! checkCelPrepared(caller) ) return 0 ;
  if (componentIndex==COMPONENT_UNDEFINED)
   {
    std::cerr<<caller<<" undefined component"<<std::endl ;
    return 0 ;
   }
  CelEventComponent * comp = getComponent(componentIndex) ;
  UShort_t treeIndex = comp->currentEntryIndex().treeIndex() ;
  return comp->currentFileSet().getFileName(treeIndex) ;
 }

const TObjString * CompositeEventList::treeName( UInt_t componentIndex ) const
 {
  TString caller ("[CompositeEventList::fileName]") ;
  if ( ! checkCelOk(caller) ) return 0 ;
  if ( ! checkCelPrepared(caller) ) return 0 ;
  if (componentIndex==COMPONENT_UNDEFINED)
   {
    std::cerr<<caller<<" undefined component"<<std::endl ;
    return 0 ;
   }
  CelEventComponent * comp = getComponent(componentIndex) ;
  UShort_t treeIndex = comp->currentEntryIndex().treeIndex() ;
  return comp->currentFileSet().getTreeName(treeIndex) ;
 }

void CompositeEventList::setDataAddress
 ( const TString & componentName,
   const TString & branchName,
   void * address )
 {
  TString caller ("CompositeEventList::setDataAddress") ;
  if ( ! checkCelPrepared(caller) ) return ;
    
  UInt_t index = componentIndex(componentName) ;
  if (index==COMPONENT_UNDEFINED)
   {
    std::cerr<<"[CompositeEventList::setDataAddress] undefined component"<<std::endl ;
    return ;
   }
  CelEventComponent * comp = getComponent(index) ;
  if (comp==0)
   {
    std::cerr<<"[CompositeEventList::setDataAddress] component not found"<<std::endl ;
    return ;
   }
  comp->setDataAddress(branchName,address) ;
 }

// If successful, returns the number of bytes read
// Failure codes:
//   -1 to -4 -> see shallowRead()
//   -5 -> Can't find a component
//   -6 -> Can't read a component input tree
Int_t CompositeEventList::deepRead( Long64_t eventIndex )
 {
  TString caller ("CompositeEventList::deepRead") ;
  if ( ! checkCelOk(caller) ) return -1 ;
  if ( ! checkCelPrepared(caller) ) return -1 ;
    
  // if there is a change of set during the shallow read
  // we must reset all the caches
  Long64_t oldSetIndex = _currentLink.setIndex() ;
  Int_t total = shallowRead(eventIndex) ;
  if (oldSetIndex!=_currentLink.setIndex())
   {
    std::vector<CelEventComponent *>::iterator itr ;
    for ( itr = _compList.begin() ; itr != _compList.end() ; itr++ )
     { (*itr)->resetSet() ; }
   }
  
  if ( total < 0 ) return total ;
  std::vector<CelEventComponent *>::iterator itr ;
  for ( itr = _compList.begin() ; itr != _compList.end() ; itr++ )
   {
    CelEventComponent * comp = *itr ;
    if ( comp == 0 ) return -5 ;
    Int_t nRead = comp->deepRead() ;
    if ( nRead < 0 ) return -6 ;
//    TTree * t = comp->getTree() ;
//    itr->first = t ;
    total += nRead ;
   }
  return total ;
 }

// Build event and data chains
TChain * CompositeEventList::newChains( TObjArray * chainList, Bool_t setFriends )
 {
  TString caller ("CompositeEventList::newChains") ;
  if ( ! checkCelOk(caller) ) return 0 ;
  if ( ! checkCelPrepared(caller) ) return 0 ;
    
  // build main chain
  TFile * f = rootUtil::getFile(*_entryTree) ;
  if (f==0) return 0 ;    
  TChain * mainChain = new TChain("EventLinks") ;  
  Int_t check = mainChain->Add(f->GetName()) ;
  if ( check < 0 )
   {
    delete mainChain ;
    return 0 ;
   }

  // build components chains
  UInt_t i ;
  for ( i=0 ; i<numComponents() ; i++ )
   {
    std::cout<<"Building chain for "<<_compNames[i]<<std::endl ;
    TChain * c = newChain(i) ;
    if (c==0) { return 0 ; }
    c->SetBranchStatus("*",0) ; // ???
    TVirtualIndex * vIdx = new CelIndex(*this,_compNames[i],c) ;
    if (vIdx==0) { return 0 ; }     
    c->SetTreeIndex(vIdx) ;
    if (chainList!=0) { chainList->AddLast(c) ; }
    if ( setFriends )
     { mainChain->AddFriend(c,_compNames[i].Data()) ; }
  }
  
  // end
  return mainChain ;
 }




//====================================================================
// NOT OK... merge use-case
//====================================================================


CompositeEventList::CompositeEventList
 ( TTree * entryTree, TTree * linkTree, TTree * fileTree )
 : _currentFile(0), _entryTree(entryTree), _linkTree(linkTree), _fileTree(fileTree)
 {
  // prepare components
  discoverComponents() ;
  
  // prepare and connect new offset tree
  _offsetTree = new TTree(__offsetTreeName.c_str(),__offsetTreeName.c_str()) ;
  Int_t check = makeCelBranches(0,0,0,_offsetTree,32000) ;
  if ( check < 0 )
   {
    // Failed to make branches for some reason.  
    // Already warned.  clean up and return NULL
    deleteCurrentFile() ;
    return ;
   }
  
  // attach all trees, and insert first 0 offsets
  attachToTree(_entryTree,_linkTree,_fileTree,_offsetTree) ;
  std::vector< CelEventComponent * >::iterator compItr ;
  for ( compItr = _compList.begin() ; compItr != _compList.end() ; compItr++ )
   {
    CelEventComponent * comp = *compItr ;
    if ( comp == 0 )
     { std::cerr<<"[CompositeEventList::CompositeEventList(,,,)] problem with components"<<std::endl ; return ; }
    comp->currentFileOffset().reset() ;
   }
  _offsetTree->Fill() ;
  
  // compute the offsets
  _fileTree->SetBranchStatus("*",kFALSE) ;
  _fileTree->SetBranchStatus("*_Set_Entries",kTRUE) ;
  Long64_t nEntries, iSet, nSet = _fileTree->GetEntries() ;
  for ( iSet = 0 ; iSet < (nSet-1) ; iSet++ )
   {
    Int_t nB = _fileTree->GetEntry(iSet) ;
    if ( nB < 0 )
     { std::cerr << "Failed to read entry in file tree " << std::endl ; }
    for ( compItr = _compList.begin() ; compItr != _compList.end() ; compItr++ )
     {
      CelEventComponent * comp = *compItr ;
      if ( comp == 0 )
       { std::cerr<<"[CompositeEventList::CompositeEventList(,,,)] problem with components"<<std::endl ; return ; }
      nEntries = comp->currentFileSet().entries() ;
      comp->currentFileOffset().increment(nEntries) ;
     }
    _offsetTree->Fill() ;
   }
 }



//====================================================================
// print and debug
//====================================================================


// Dump a set of event component pointers and the list of TTree where they live
void CompositeEventList::printInfo
 ( const char * options, UInt_t nEvent, UInt_t startEvent )
 {
  if (nEvent==0)
   { nEvent = numEvents() ; }
  UInt_t lastEvt = startEvent+nEvent ;
        
  TString caller ("CompositeEventList::printInfo") ;
  if ( ! checkCelPrepared(caller) ) return ;
  
  // init widths
  _treeNumberWidth = 0 ;
  _entryNumberWidth = 0 ;
  UInt_t nEvents = lastEvt ;
  _eventNumberWidth = 1 ; // width
  while (nEvents>10)
   { ++_eventNumberWidth ; nEvents /= 10 ; }
  UInt_t nSets = numFileAndTreeSets() ;
  _setNumberWidth = 1 ;
  while (nSets>10)
   { ++_setNumberWidth ; nSets /= 10 ; }
      
  // Components
  std::cout << "Printing Component Names:" ;
  _componentNameWidth = 0 ;
  TString componentName ;
  std::vector<CelEventComponent*>::const_iterator itr ;
  for ( itr = _compList.begin() ; itr != _compList.end() ; itr++ )
   {
    const CelEventComponent * comp = *itr ;
    if ( comp == 0 ) throw "?!?" ;
    componentName = comp->componentName() ;
    std::cout << " " << componentName ;
    if (((unsigned int)componentName.Length())>_componentNameWidth)
     { _componentNameWidth = componentName.Length() ; }
   } 
  std::cout << std::endl ;
        
  // Events
  std::cout << "Printing Events Info: " << std::endl ;
  UInt_t iEvt ;
  for ( iEvt = startEvent ; iEvt < lastEvt ; iEvt++ )
   {
    shallowRead(iEvt) ;
    printEventInfo(options) ;
   }
  
  // Sets
  std::cout << "Printing Sets Info: " << std::endl ;
  Long64_t iSet = -1 ;
  for ( iEvt = startEvent ; iEvt < lastEvt ; iEvt++ )
   {
    shallowRead(iEvt) ;
    if (currentSetIndex()!=iSet)
     {
      iSet = currentSetIndex() ;
      printSetInfo(options) ;
     }
   }
 }

// Dump a single event on one line
void CompositeEventList::printEventInfo( const char * options )
 {
  // should be done once for all somewhere else...
  Long64_t nEvents = numEvents() ;
  Long64_t nSets = numFileAndTreeSets() ;
  unsigned int wEvents = 1 ; // width
  while (nEvents>10) { ++wEvents ; nEvents /= 10 ; }
  unsigned int wSets = 1 ;
  while (nSets>10) { ++wSets ; nSets /= 10 ; }
  
  // print info
  std::cout
    << "Event " << setw(_eventNumberWidth) << std::right << _currentLink.eventIndex() << " : "
    << "Set " << setw(_setNumberWidth) << std::right << _currentLink.setIndex() ;
  std::vector<CelEventComponent*>::const_iterator itr ;
  for ( itr = _compList.begin() ; itr != _compList.end() ; itr++ )
   {
    const CelEventComponent * comp = *itr ;
    if ( comp == 0 ) return ; // ??
    std::cout << ", ";
    std::cout << comp->componentName() << ' ';
    comp->printEventInfo(options) ;
   } 
  std::cout << std::endl ;
 }

// Dump the list of TTree where our events live
void CompositeEventList::printSetInfo( const char * options )
 {
  TString setPrefix("Set ") ;
  TString name ;
  TString prefix ;
  setPrefix += currentSetIndex() ;
  std::vector<CelEventComponent*>::const_iterator itr ;
  for ( itr = _compList.begin() ; itr != _compList.end() ; itr++ )
   {
    const CelEventComponent * comp = *itr ;
    if ( comp == 0 ) return ; // ??
    name = comp->componentName() ;
    while (name.Length()<_componentNameWidth)
     { name += ' ' ; }
    prefix = setPrefix + ", " + name ;
    comp->printSetInfo(options,prefix) ;
  }   
}




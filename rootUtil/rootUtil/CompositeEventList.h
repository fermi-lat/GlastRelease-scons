//
//========================================================
// File: CompositeEventList.h
// Authors:
//   DC, David Chamont, LLR, chamont@llr.in2p3.fr
//   EC, Eric Charles , SLAC, echarles@slac.stanford.edu
//========================================================
//



#ifndef CompositeEventList_h
#define CompositeEventList_h



#include "CelEventLink.h"
class CelEventComponent ;
class CelIndex ;

#include <TObject.h>
#include <TString.h>

class TTree ;
class TFile ;
class TChain ;
class TVirtualIndex ;
class TCollection ;
class TObjArray ;
class TObjString ;
class RuChain ;

#include <map>
#include <vector>

//
// The class CompositeEventList permits to represent in memory
// a list of composite events, where each event is described as
// a set of entries in several separate data trees. Each data tree
// contains a different kind of data, and is also called a "component".
//
// An instance of CompositeEventList can be read from and/or
// written to a ROOT file. We will also implement later the construction
// of a CompositeEventList on-the-fly, so to support the good old
// way to read/write events as was done by former RootIo writers
// and readers.
//
// The information is stored within 5 internal trees
// 
//   TO BE DONE : EventIDs: 3 branches, 1 entry per event
//      Event_ProductionID (char *)
//      Event_RunID        (Long64_t)
//      Event_EventID      (Long64_t)
// 
//   EventEntries: 2 Branches per component, 1 entry per event
//      <Comp>_Event_EntryIndex (Long64_t) : Index of the entry in the tree it lives on
//      <Comp>_Event_TreeIndex  (Short_t)  : Index of the tree in the associated set of files and trees
//
//   EventLinks: 2 branches, 1 entry per event
//      Event_SetIndex  (Long64_t) : Index of the associated set of files and trees
//      Event_SetOffset (Long64_t) : Number of events in previous sets
// 
//   FileAndTreeSets: 5 Branches per component, 1 entry per set
//      <Comp>_Set_Size        -> Number of Files and Trees in the current set
//      <Comp>_Set_FileNames   -> TClonesArray<TObjString> with the file names
//      <Comp>_Set_TreeNames   -> TClonesArray<TObjString> with the tree names (should all be the same ?)
//      <Comp>_Set_TreeOffsets -> TArrayL64 with the offsets for each Tree in the current set
//      <Comp>_Set_Entries     -> Total number of entries in all trees from the current set
//
// Note : the components can be given when creating the CEL objet,
// and/or afterwards with declareComponent() method.
//
// Note : the pair (RunID/EventID) is expected to be unique for GLAST real data. For
// what concerns the simulation data, the triplet (ProductionID/RunID/EventID) should
// be unique.
//
// CAUTION : the fillEvent() method rely on TTree::GetReadEntry() to know
// which is the current entry in each TTree. If you should just filled it,
// TTree::GetReadEntry() could not return the last filled entry index, and
// you should call tree->TTree::LoadTree(tree->GetEntries()-1) on each tree
// before calling CompositeEventList::fillEvent().

class CompositeEventList : public TObject
 {

  public :
   
    // Construction and components
    explicit CompositeEventList
     ( const TString & celFileName = "",
       const TString & options = "READ",
       const TObjArray * componentNames = 0 ) ;
    UInt_t declareComponent( const TString & name ) ;
    Bool_t isOk() ;
    ~CompositeEventList() ;

    // OK! Write
    Bool_t fillEntry(  const TString & componentName, TTree * ) ;  // components one by one
    Bool_t fillEntry(  UInt_t componentIndex, TTree * ) ;  // components one by one
    Long64_t fillEvent() ;                                // components one by one
    Long64_t fillEvent( const TObjArray & ) ;            // all components together
    Long64_t fillEvent( const std::vector<TTree*> & ) ;  // all components together
    Long64_t fillEvent( const std::vector<TChain*> & ) ; // all components together
    Long64_t fillFileAndTreeSet() ;
    void writeAndClose() ;
    
    // OK! Shallow read (do not handle the user data trees)
    static const UInt_t COMPONENT_UNDEFINED ;
    UInt_t componentIndex( const TString &  componentName ) const ;
    RuChain * newChain( UInt_t componentIndex ) const ;
    Int_t shallowRead( Long64_t eventIndex ) ;
    Long64_t entryIndex( UInt_t componentIndex ) const ;
    const TObjString * fileName( UInt_t componentIndex ) const ;
    const TObjString * treeName( UInt_t componentIndex ) const ;
    
    // OK! Deep read (handle the user data trees)
    void setDataAddress
     ( const TString & componentName,
       const TString & branchName,
       void * address ) ;
    Int_t deepRead( Long64_t eventIndex ) ;   
    
    // NOT OK! Read through friend chains
    TChain * newChains( TObjArray * chainList = 0, Bool_t setFriends = kTRUE ) ;
    
    // Various Accessors
    const TString & fileName() const ;
    UInt_t numComponents() const ;
    const TString & componentName( UInt_t componentIndex ) const ;
    Long64_t numEvents() const ;
    Long64_t currentEventIndex() const ;
    Long64_t numFileAndTreeSets() const ;
    Long64_t currentSetIndex() const ;
    
    /// PRINTING
    // Dump a set of event component pointers and the list of TTree where they live
    void printInfo( const char * options ="", UInt_t nEvent=0, UInt_t startEvent=0) ;
    // Dump a single event
    void printEventInfo( const char * options ="" ) ;
    // Dump the list of TTree where our events live
    void printSetInfo( const char * options ) ;

    /// FOR CelUtil
    TTree * entryTree() { return _entryTree ; }
    TTree * linkTree() { return _linkTree ; }
    TTree * fileTree() { return _fileTree ; }
    TTree * offsetTree() { return _offsetTree ; }
    CompositeEventList( TTree * eventTree, TTree * fileTree, TTree * entryTree ) ;

    /// FOR CelIndex
    CelEventComponent * getComponent( const TString & name ) const ;
  
     
  private :
	  
//    typedef std::pair<TTree*,CelEventComponent*> CelTreeAndComponent ;
//	  typedef std::vector<CelTreeAndComponent> CelTreesAndComponents ;

    // Disable copying and assignment
    CompositeEventList( const CompositeEventList & ) ;
    CompositeEventList & operator=( const CompositeEventList & ) ;

    // Manipulation of cel internal trees
    void deleteCurrentFile() ;
    Bool_t checkCelOk( const TString & caller ="CompositeEventList::?" ) const ;
    Bool_t checkCelPrepared( const TString & caller ="CompositeEventList::?" ) const ;
    Bool_t checkCelPrepared( const TString & caller ="CompositeEventList::?" ) ;
    Bool_t checkCelTrees( const TString & caller ="CompositeEventList::?" ) ;
    Bool_t checkCelTree( const TString & caller, TTree *, const TString & name, Bool_t error ) ;
    Int_t makeCelBranches( TTree * entryTree, TTree * linkTree, TTree * fileTree, TTree * offsetTree, Int_t bufsize = 32000) const;
    Int_t attachToTree( TTree * entryTree, TTree * linkTree, TTree * fileTree, TTree * offsetTree ) const ;

    // Components utilities
    UInt_t CompositeEventList::declareComponents( const TObjArray * componentNames ) ;
    UInt_t discoverComponents() ;
    CelEventComponent * getComponent( UInt_t index ) const ;
    // UNUSED ? Get a Tree that is being read by a deepRead
    TTree * getTree( UInt_t componentIndex ) const ;
    TTree * getTree( const TString & componentName ) const ;

    // specific construction methods
    void prepareRead() const ;
    void prepareRecreate() ;
    
    // cel data
    mutable Bool_t  _isOk ;
    mutable Bool_t  _isPrepared ;
    TString _fileName ;
    TString _openingOptions ;
    TFile * _currentFile ;
    TTree * _entryTree ; 
    TTree * _linkTree ; 
    TTree * _fileTree ; 
    TTree * _offsetTree ; 
    
    // current event info
    CelEventLink _currentLink ;

    // components stuff
    std::vector<CelEventComponent*> _compList ;  //! vector of components
    std::vector<TString> _compNames ;   //! names of components
    std::map<TString,UInt_t> _compMap ; //! components lookup map to get the index in the list above from the name

    // internal print stuff
    Ssiz_t _componentNameWidth ;
    unsigned int _setNumberWidth ;
    unsigned int _eventNumberWidth ;
    unsigned int _treeNumberWidth ;
    unsigned int _entryNumberWidth ;
    
    ClassDef(CompositeEventList,0)

 } ;


#endif



#ifndef TestData_h
#define TestData_h


#include "ComponentsInfo.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#  include <TObject.h>
#  include <TString.h>
#  include <TFile.h>
#  include <TChain.h>
#  include <TObjArray.h>
#  include <TRandom.h>
#  include <Riosfwd.h>
#  include <Riostream.h>
#  include <vector>
#endif


//==================================================
// Dummy labels
//==================================================


class TestDigiLabel {} ;
class TestReconLabel {} ;



//==================================================
// Dummy components
//==================================================


class TestAbstractComponent
 {
  public :
    virtual const ComponentInfo * getInfo() const =0 ;
  protected :
    static const ComponentsInfo * allInfo() ;
 } ;
        
template <class Label>
class TestComponent : public TestAbstractComponent
 {
  public :
    static const ComponentInfo * info() ;
    virtual const ComponentInfo * getInfo() const
     { return info() ; }	
 } ;
 
typedef TestComponent<TestDigiLabel> TestDigiComponent ;
typedef TestComponent<TestReconLabel> TestReconComponent ;


  
//==================================================
// Dummy data classes
//==================================================

typedef UInt_t IdType ;

class TestAbstractData
 {
  public :
    TestAbstractData( IdType runID =0, IdType eventID =0, Bool_t isValid =kTRUE )
     : runID_(runID), eventID_(eventID), isValid_(isValid) {}
    IdType getRunID() const { return runID_ ; }
    IdType getEventID() const { return eventID_ ; }
    Bool_t isValid() const { return isValid_ ; }
    void setRunID( IdType runID ) { runID_ = runID ; }
    void setEventID( IdType eventID ) { eventID_ = eventID ; }
    void setIsValid( Bool_t isValid ) { isValid_ = isValid ; }
    void incrementRunID() { ++runID_ ; }
    void incrementEventID() { ++eventID_ ; }
    virtual const TString & getComponentName() const =0 ;
  private :
    IdType runID_ ;
    IdType eventID_ ;
    Bool_t isValid_ ;
    ClassDef(TestAbstractData,1)
 } ;
 
std::ostream & operator<<( std::ostream &, const TestAbstractData & ) ;

template <class Label>
class TestData : public TestAbstractData
 {
  public :
    TestData( IdType runID =0, IdType eventID =0 )
     : TestAbstractData(runID,eventID) {}
    virtual const TString & getComponentName() const
     { return TestComponent<Label>::info()->componentName ; }
  private :
    ClassDef(TestData,1)
 } ;

typedef TestData<TestDigiLabel> TestDigi ;
typedef TestData<TestReconLabel> TestRecon ;



//==================================================================
// Simple Writer : each component and run is generated independantly
//==================================================================


template <class Label>
int testWrite( char * baseName, IdType runId, IdType firstEvent, IdType lastEvent, TRandom * random, Bool_t makeIndex = kFALSE )
 {   
  Int_t buffer = 64000 ;
  Int_t splitLevel = 99 ;
  
  TString fileName = baseName ;
  fileName += TestComponent<Label>::info()->componentName ;
  fileName += "." ;
  fileName += runId ;
  fileName += "." ;
  if (firstEvent<10)
   { fileName += "0" ; }
  fileName += firstEvent ;
  fileName += "-" ;
  if (lastEvent<10)
   { fileName += "0" ; }
  fileName += lastEvent ;
  fileName += ".root" ;

  std::cout
    <<"[testWrite] Creating file "<<fileName
    <<std::endl ;
  
  TFile * f = new TFile(fileName,"RECREATE") ;
  TTree * t = new TTree(TestComponent<Label>::info()->treeName,TestComponent<Label>::info()->treeName) ;
    
  std::cout<<"[testWrite] Creating branches"<<std::endl ;
  TestData<Label> * componentEntry = new TestData<Label>(runId,firstEvent) ;
  t->Branch(TestComponent<Label>::info()->topBranchName,TestComponent<Label>::info()->topBranchType,&componentEntry,buffer,splitLevel) ;
  
  std::cout
    <<"[testWrite] Filling "
    <<TestComponent<Label>::info()->treeName<<" tree" << std::endl ;
  
  IdType ievent = firstEvent, eventID  ;
  while ( ievent <= lastEvent )
   {
    if (random&&random->Uniform()>0.5)
     {
      eventID = componentEntry->getEventID() ;
      componentEntry->setIsValid(kFALSE) ;
      t->Fill() ;
      componentEntry->setIsValid(kTRUE) ;
     }
    else
     {
      t->Fill();
      componentEntry->incrementEventID() ;
      ievent++ ;
     }
   }
  std::cout << "[testWrite] File filled with " << (lastEvent-firstEvent+1) << " events" << std::endl ;
  
  if (makeIndex==kTRUE)
   {
    t->BuildIndex
     ( TestComponent<Label>::info()->runIdBranchName,
       TestComponent<Label>::info()->eventIdBranchName ) ;
    std::cout << "[testWrite] Index made" << std::endl ;
   }

  std::cout << "[testWrite] Writing and closing file" << std::endl ;
  delete componentEntry ;
  f->Write() ;
  f->Close() ;
  delete f ;
  return 0 ;
 }



//============================================================
// Simple Reader : we want to read several components together
//============================================================


class TestReader
 {
  public :
	  
    TestReader() {}
    void add( const char * baseName, const ComponentInfo * component ) ;
    
    void showByComponent() const ;
    void showByEvent() ;
    
    void resetEvent() ;
    bool nextEvent() ;
    void showEvent() const ;
    
    typedef std::vector<TChain*> RootForest ;
    const RootForest & getForest() const
     { return forest_ ; }

  private :
	  
    TObjArray readers_ ;
    RootForest forest_ ;
    IdType currentEvent_ ;
	
    struct BranchReader : public TObject
     {
      const ComponentInfo * component_ ;
      TChain * chain_ ;
      Long64_t nextEntry_ ;
      TestAbstractData * data_ ;
     } ;

 } ;



//==================================================
// End
//==================================================


#endif



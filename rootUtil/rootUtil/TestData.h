

#ifndef TestData_h
#define TestData_h


#include <TObject.h>
#include <TString.h>
#include <TFile.h>
#include <TChain.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <Riosfwd.h>



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
    virtual const TString & getName() =0 ;
    virtual const TString & getTreeName() =0 ;
    virtual const TString & getBranchName() =0 ;
    virtual const TString & getDataTypeName() =0 ;
    //virtual TestAbstractData * newData( Long64_t runID =0, Long64_t eventID =0 ) =0 ;
 } ;
	
template <class Label>
class TestComponent : public TestAbstractComponent
 {
  public :
	  
	static  const TString & name() { return componentName__ ; }
	static  const TString & treeName() { return treeName__ ; }
	static  const TString & branchName() { return branchName__ ; }
	static  const TString & dataTypeName() { return dataTypeName__ ; }
	
	virtual const TString & getName() { return componentName__ ; }
	virtual const TString & getTreeName() { return treeName__ ; }
	virtual const TString & getBranchName() { return branchName__ ; }
	virtual const TString & getDataTypeName() { return dataTypeName__ ; }
	
	//virtual TestAbstractData * newData( Long64_t runID =0, Long64_t eventID =0 )
	// { return new TestData<Label>(runID,eventID) ; }
	
  private :
	  
	static TString componentName__ ;
	static TString treeName__ ;
	static TString branchName__ ;
	static TString dataTypeName__ ;
 } ;
 
typedef TestComponent<TestDigiLabel> TestDigiComponent ;
typedef TestComponent<TestReconLabel> TestReconComponent ;


  
//==================================================
// Dummy data classes
//==================================================


class TestAbstractData
 {
  public :
    TestAbstractData( Long64_t runID =0, Long64_t eventID =0)
     : runID_(runID), eventID_(eventID) {}
	Long64_t getRunID() const { return runID_ ; }
	Long64_t getEventID() const { return eventID_ ; }
	void setRunID( Long64_t runID ) { runID_ = runID ; }
	void setEventID( Long64_t eventID ) { eventID_ = eventID ; }
	void incrementRunID() { ++runID_ ; }
	void incrementEventID() { ++eventID_ ; }
	virtual const TString & getComponentName() const =0 ;
  private :
	Long64_t runID_ ;
	Long64_t eventID_ ;
    ClassDef(TestAbstractData,1)
 } ;
 
std::ostream & operator<<( std::ostream &, const TestAbstractData & ) ;

template <class Label>
class TestData : public TestAbstractData
 {
  public :
    TestData( Long64_t runID =0, Long64_t eventID =0 )
     : TestAbstractData(runID,eventID) {}
	virtual const TString & getComponentName() const
	 { return TestComponent<Label>::name() ; }
  private :
    ClassDef(TestData,1)
 } ;

typedef TestData<TestDigiLabel> TestDigi ;
typedef TestData<TestReconLabel> TestRecon ;



//==================================================================
// Simple Writer : each component and run is generated independantly
//==================================================================


template <class Label>
int testWrite( char * baseName, Long64_t runId, Long64_t firstEvent, Long64_t lastEvent, TRandom & random )
 {   
  Int_t buffer = 64000 ;
  Int_t splitLevel = 1 ;
  
  TString fileName = baseName ;
  fileName += TestComponent<Label>::name() ;
  fileName += "." ;
  fileName += runId ;
  fileName += "." ;
  fileName += firstEvent ;
  fileName += "-" ;
  fileName += lastEvent ;
  fileName += ".root" ;

  std::cout
    <<"[testWrite] Creating file "<<fileName
    <<std::endl ;
  
  TFile * f = new TFile(fileName,"RECREATE") ;
  TTree * t = new TTree(TestComponent<Label>::treeName(),TestComponent<Label>::treeName()) ;
    
  std::cout << "[testWrite] Creating branches" << std::endl ;
  TestData<Label> * componentEntry = new TestData<Label>(runId,firstEvent) ;
  t->Branch(TestComponent<Label>::branchName(),TestComponent<Label>::dataTypeName(),&componentEntry,buffer,splitLevel) ;
  
  std::cout
    <<"[testWrite] Filling "
    <<TestComponent<Label>::treeName()<<" tree" << std::endl ;
  
  Long64_t ievent = firstEvent, eventID  ;
  while ( ievent <= lastEvent )
   {
	if (random.Uniform()>0.5)
	 {
	  eventID = componentEntry->getEventID() ;
	  componentEntry->setRunID(-1) ;
	  componentEntry->setEventID(-1) ;
      t->Fill() ;
	  componentEntry->setEventID(eventID) ;
	  componentEntry->setRunID(runId) ;
	 }
	else
	 {
      t->Fill();
      componentEntry->incrementEventID() ;
      ievent++ ;
	 }
   }
  std::cout << "[testWrite] file filled with " << (lastEvent-firstEvent+1) << " events" << std::endl ;
  delete componentEntry ;

  std::cout << "[testWrite] Writing and closing file" << std::endl ;
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
    void add( const char * baseName, TestAbstractComponent * component ) ;
    
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
	Long64_t currentEvent_ ;
	
	
	struct BranchReader : public TObject
     {
      TestAbstractComponent * component_ ;
      TChain * chain_ ;
      Long64_t nextEntry_ ;
      TestAbstractData * data_ ;
     } ;

 } ;



//==================================================
// End
//==================================================


#endif



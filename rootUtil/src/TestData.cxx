

#include <rootUtil/TestData.h>
#include <Riostream.h>



//==================================================
// Dummy components
//==================================================


// Below, one cannot use the typedefs for the dataTypeName__
// values, because it is reused for TTree:Branch(), which
// does not know about typedefs.

TString TestDigiComponent::componentName__ = "Digi" ;
TString TestDigiComponent::treeName__ = "TestDigiTree" ;
TString TestDigiComponent::branchName__ = "TestDigiBranch" ;
TString TestDigiComponent::dataTypeName__ = "TestData<TestDigiLabel>" ;
ClassImp(TestDigi)

TString TestReconComponent::componentName__ = "Recon" ;
TString TestReconComponent::treeName__ = "TestReconTree" ;
TString TestReconComponent::branchName__ = "TestReconBranch" ;
TString TestReconComponent::dataTypeName__ = "TestData<TestReconLabel>" ;
ClassImp(TestRecon)


  
//==================================================
// Dummy data classes
//==================================================


std::ostream & operator<<( std::ostream & os, const TestAbstractData & data )
 { return (os<<data.getComponentName()<<" "<<data.getRunID()<<"|"<<data.getEventID()) ; }



//============================================================
// Simple Reader : we want to read several components together
//============================================================


void TestReader::add( const char * baseName, TestAbstractComponent * component )
 {
  BranchReader * reader = new BranchReader ;
  
  reader->component_ = component ;
  reader->chain_ = new TChain(component->getTreeName()) ;
  
  TString fileNames = baseName ;
  fileNames += component->getName() ;
  fileNames += "*.root" ;
  reader->chain_->Add(fileNames) ;
  
  reader->data_ = 0 ;
  reader->chain_->SetBranchAddress
   (component->getBranchName(),&(reader->data_)) ;
  
//  if (readers_.GetEntries()==0)
//   { nbEvents_ = reader->chain_->GetEntries() ; }
//  else if (nbEvents_!=reader->chain_->GetEntries())
//   { std::cout<<"[TestReader::add] inconsistent chain sizes"<<std::endl ; }

  readers_.Add(reader) ;
 }
	
void TestReader::showByComponent()
 {
  std::cout<<"[TestReader] start looping on components"<<std::endl ;
  TIter next(&readers_) ;
  BranchReader * reader = 0 ;
  while ((reader=(BranchReader *)next()))
   {
    std::cout<<"[TestReader] Component "<<reader->component_->getName() ;
    TString separator(" : ") ;
    UInt_t nbEvents = reader->chain_->GetEntries() ;
    UInt_t ievent ;
    for ( ievent=0 ; ievent<nbEvents ; ++ievent )
     {
      reader->chain_->GetEntry(ievent) ;
      std::cout
        <<separator<<reader->data_->getRunID()
        <<'|'<<reader->data_->getEventID() ;
      separator = ", " ;
     }
    std::cout<<std::endl ;
   }
 }

void TestReader::showByEvent()
 {
  // reset entries
  TIter next(&readers_) ;
  BranchReader * reader = 0 ;
  while ((reader=(BranchReader *)next()))
   { reader->nextEntry_ = 0 ; }

  // loop, ignoring -1 values
  std::cout<<"[TestReader] start looping on events"<<std::endl ;
  bool finished = false ;
  UInt_t ievent = 0 ;
  while ( ! finished )
   {
    std::cout<<"[TestReader] Event "<<ievent ;
    TString separator(" : ") ;
    TIter next(&readers_) ;
    BranchReader * reader = 0 ;
    while ((reader=(BranchReader *)next()))
     {
      while (reader->nextEntry_<reader->chain_->GetEntries())
       {
        reader->chain_->GetEntry(reader->nextEntry_++) ;
        if (reader->data_->getRunID()!=-1)
         { break ; }
       }
      std::cout<<separator<<*(reader->data_) ;
      separator = ", " ;
      if (reader->nextEntry_==reader->chain_->GetEntries())
       { finished = true ; }
     }
    ++ievent ;
    std::cout<<std::endl ;
   }
 }


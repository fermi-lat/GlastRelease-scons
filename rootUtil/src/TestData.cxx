

#include <rootUtil/TestData.h>
#include <Riostream.h>
#include <cassert>



//==================================================
// Dummy components
//==================================================


// Below, one cannot use the typedefs for the dataTypeName__
// values, because it is reused for TTree:Branch(), which
// does not know about typedefs.

template <> const TString & TestDigiComponent::name()
 { static TString text = "Digi" ; return text ; }
template <> const TString & TestDigiComponent::treeName()
 { static TString text = "TestDigiTree" ; return text ; }
template <> const TString & TestDigiComponent::branchName()
 { static TString text = "TestDigiBranch" ; return text ; }
template <> const TString & TestDigiComponent::dataTypeName()
 { static TString text = "TestData<TestDigiLabel>" ; return text ; }
ClassImp(TestDigi)

template <> const TString & TestReconComponent::name()
 { static TString text = "Recon" ; return text ; }
template <> const TString & TestReconComponent::treeName()
 { static TString text = "TestReconTree" ; return text ; }
template <> const TString & TestReconComponent::branchName()
 { static TString text = "TestReconBranch" ; return text ; }
template <> const TString & TestReconComponent::dataTypeName()
 { static TString text = "TestData<TestReconLabel>" ; return text ; }
ClassImp(TestRecon)


  
//==================================================
// Dummy data classes
//==================================================


std::ostream & operator<<( std::ostream & os, const TestAbstractData & data )
 {
  return
   (
    os
      <<data.getComponentName()<<" "
      <<setw(2)<<std::right<<data.getRunID()<<"|"
      <<setw(2)<<std::left<<data.getEventID()
   ) ;
 }



//==========================================
// Simple Reader : déclaration of components
//==========================================


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
  forest_.push_back(reader->chain_) ;
 }
	


//==============================================================
// Simple Reader : show what is within each separate component
//==============================================================


void TestReader::showByComponent() const
 {
  std::cout<<"[TestReader] Looping on components" ;
  TIter next(&readers_) ;
  BranchReader * reader = 0 ;
  while ((reader=(BranchReader *)next()))
   {
    TString prefix("\n[TestReader] Component ") ;
    prefix += reader->component_->getName() ;
    prefix += " : " ;
    TString separator ;
    
    UInt_t nbEvents = reader->chain_->GetEntries() ;
    UInt_t ievent ;
    for ( ievent=0 ; ievent<nbEvents ; ++ievent )
     {
      if (ievent%8)
       { separator = ", " ; }
      else
       { separator = prefix ; }
      reader->chain_->GetEntry(ievent) ;
      if ((reader->data_->getRunID()==-1)&&(reader->data_->getEventID()==-1))
       { std::cout<<separator<<" *|* " ; }
      else
       {
        std::cout
          <<separator<<setw(2)<<std::right<<reader->data_->getRunID()
          <<'|'<<setw(2)<<std::left<<reader->data_->getEventID() ;
       }
     }
   }
  std::cout<<std::endl ;
 }



//==========================================
// Simple Reader : access events one by one
//==========================================


void TestReader::resetEvent()
 {
  TIter next(&readers_) ;
  BranchReader * reader = 0 ;
  while ((reader=(BranchReader *)next()))
   { reader->nextEntry_ = 0 ; }
  currentEvent_ = -1 ;
 }

bool TestReader::nextEvent()
 {
  // search next component entries, ignoring -1 values
  bool result = true ;
  TIter next(&readers_) ;
  BranchReader * reader = 0 ;
  while ((reader=(BranchReader *)next()))
   {
    if (reader->nextEntry_==reader->chain_->GetEntries())
     { result = false ; }
    else
     {
      while (reader->nextEntry_<reader->chain_->GetEntries())
       {
        reader->chain_->GetEntry(reader->nextEntry_++) ;
        if ((reader->data_->getRunID()!=-1)&&
    	    (reader->data_->getEventID()!=-1))
         { break ; }
       }
      if ((reader->data_->getRunID()==-1)||(reader->data_->getEventID()==-1))
       { result = false ; }
     }
   }
  if ( result == true )
   { currentEvent_++ ; }
  else
   { resetEvent() ; }
  return result ;
 }

void TestReader::showEvent() const
 {
  assert(currentEvent_!=-1) ; 
  
  std::cout<<"[TestReader] Event "<<setw(2)<<std::right<<currentEvent_ ;
  TString separator(" : ") ;
    
  TIter next(&readers_) ;
  BranchReader * reader = 0 ;
  while ((reader=(BranchReader *)next()))
   {
    std::cout<<separator<<*(reader->data_) ;
    separator = ", " ;
   }
  std::cout<<std::endl ;
 }

void TestReader::showByEvent()
 {
  resetEvent() ;
  std::cout<<"[TestReader] Looping on events"<<std::endl ;
  while ( nextEvent() )
   { showEvent() ; }
 }


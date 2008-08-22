

#include <rootUtil/TestData.h>
#include <rootUtil/ComponentsInfoGlast.h>

#if !defined(__CINT__) || defined(__MAKECINT__)
#  include <Riostream.h>
#endif

#include <cassert>



//==================================================
// Dummy components
//==================================================


const ComponentsInfo * TestAbstractComponent::allInfo()
 { static ComponentsInfoGlast allTheInfo ; return &allTheInfo ; }
        
template <> const ComponentInfo * TestDigiComponent::info()
 { static const ComponentInfo * myInfo = allInfo()->getInfo("tdigi") ; return myInfo ; }
//ClassImp(TestDigi)

template <> const ComponentInfo * TestReconComponent::info()
 { static const ComponentInfo * myInfo = allInfo()->getInfo("trecon") ; return myInfo ; }
//ClassImp(TestRecon)


  
//==================================================
// Dummy data classes
//==================================================


std::ostream & operator<<( std::ostream & os, const TestAbstractData & data )
 {
  return
   (
     ( data.isValid()==kTRUE ) ?
     ( os
         <<data.getComponentName()<<" "
         <<setw(2)<<std::right<<data.getRunID()<<"|"
         <<setw(2)<<std::left<<data.getEventID() ) :
     ( os<<data.getComponentName()<<" *|* " )
   ) ;
 }



//==========================================
// Simple Reader : déclaration of components
//==========================================


void TestReader::add( const char * baseName, const ComponentInfo * component )
 {
  BranchReader * reader = new BranchReader ;
  
  reader->component_ = component ;
  reader->chain_ = new TChain(component->treeName) ;
  
  TString fileNames = baseName ;
  fileNames += component->componentName ;
  fileNames += "*.root" ;
  reader->chain_->Add(fileNames) ;
  
  reader->data_ = 0 ;
  reader->chain_->SetBranchAddress
   (component->topBranchName,&(reader->data_)) ;
  
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
    prefix += reader->component_->componentName ;
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
      std::cout<<separator<<(*(reader->data_)) ;
//      
//      if (reader->data_->isValid()==kFALSE)
//       { std::cout<<separator<<" *|* " ; }
//      else
//       {
//        std::cout
//          <<separator<<setw(2)<<std::right<<reader->data_->getRunID()
//          <<'|'<<setw(2)<<std::left<<reader->data_->getEventID() ;
//       }
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
  currentEvent_ = ((IdType)-1) ;
 }

bool TestReader::nextEvent()
 {
  // search next component valid entries
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
        if (reader->data_->isValid()==kTRUE)
         { break ; }
       }
      if (reader->data_->isValid()==kFALSE)
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
  assert(currentEvent_!=((IdType)-1)) ; 
  
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


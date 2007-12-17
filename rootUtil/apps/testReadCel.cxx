
#include <rootUtil/TestData.h>
#include <rootUtil/CompositeEventList.h>
#include <rootUtil/CelManager.h>
#include "Riostream.h"

int main( int argc, char ** argv )
 {
  char * baseName = "test" ;
  int argi = 1 ;
  if (argc > 1)
   { baseName = argv[argi++] ; } 

  int sc = 0 ;
  try 
   {
	// direct reader access
    std::cout
      << "=================================\n"
      << "DIRECT READ\n"
      << std::flush ;
    TestReader reader ;
    reader.add(baseName,new TestDigiComponent) ;
    reader.add(baseName,new TestReconComponent) ;
    reader.showByComponent() ;
    reader.showByEvent() ;
    
    TString fileName = baseName ;
    fileName += ".cel.root" ;
    CompositeEventList cel(fileName) ;
    
    // internal cel print
    std::cout
      << "=================================\n"
      << "CEL Internal meta info\n"
      << std::flush ;
    cel.printInfo() ;
    
	// usual read variables
    Long64_t nEvents, iEvent ;
    TestAbstractData * digiData =0 ;
    TestAbstractData * reconData =0 ;
    TChain * digiChain =0 ;
    TChain * reconChain =0 ;

    // use of CompositeEventList::shallowRead()
    std::cout
      << "=================================\n"
      << "SHALLOW READ\n"
      << std::flush ;
    UInt_t digiIndex = cel.componentIndex(TestDigiComponent::name()) ;
    digiChain = cel.newChain(digiIndex) ;
    digiChain->SetBranchAddress(TestDigiComponent::branchName(),&digiData) ;
    UInt_t reconIndex = cel.componentIndex(TestReconComponent::name()) ;
    reconChain = cel.newChain(reconIndex) ;
    reconChain->SetBranchAddress(TestReconComponent::branchName(),&reconData) ;
    nEvents = cel.numEvents() ;
    for ( iEvent = 0 ; iEvent < nEvents ; iEvent++ )
     {
      std::cout
        <<"Event "<<setw(2)<<std::right<<iEvent
        <<" : "<<std::flush ;    	  
      cel.shallowRead(iEvent) ;
      digiChain->GetEntry(cel.entryIndex(digiIndex)) ;
      reconChain->GetEntry(cel.entryIndex(reconIndex)) ;
      if ((digiData!=0)&&(reconData!=0))
       { std::cout<<*digiData<<", "<<*reconData<<std::endl ; }
      else
       { std::cout<<"Digi and/or Recon data not read"<<std::endl ; }
     }
    delete digiChain ;
    delete reconChain ;
    digiChain = reconChain = 0 ;
    delete digiData ;
    delete reconData ;
    digiData = reconData = 0 ;
    
    
    // use of CompositeEventList::deepRead()
    std::cout
      << "=================================\n"
      << "DEEP READ\n"
      << std::flush ;
//    UInt_t digiIndex = cel.componentIndex(TestDigiComponent::name()) ;
//    digiChain = cel.newChain(digiIndex) ;
//    digiChain->SetBranchAddress(TestDigiComponent::branchName(),&digiData) ;
//    UInt_t reconIndex = cel.componentIndex(TestReconComponent::name()) ;
//    reconChain = cel.newChain(reconIndex) ;
//    reconChain->SetBranchAddress(TestReconComponent::branchName(),&reconData) ;
//    nEvents = cel.numEvents() ;
//    for ( iEvent = 0 ; iEvent < nEvents ; iEvent++ )
//     {
//      std::cout
//        <<"Event "<<setw(2)<<std::right<<iEvent
//        <<" : "<<std::flush ;    	  
//      cel.shallowRead(iEvent) ;
//      digiChain->GetEntry(cel.entryIndex(digiIndex)) ;
//      reconChain->GetEntry(cel.entryIndex(reconIndex)) ;
//      if ((digiData!=0)&&(reconData!=0))
//       { std::cout<<*digiData<<", "<<*reconData<<std::endl ; }
//      else
//       { std::cout<<"Digi and/or Recon data not read"<<std::endl ; }
//     }
//    delete digiChain ;
//    delete reconChain ;
//    digiChain = reconChain = 0 ;
//    delete digiData ;
//    delete reconData ;
//    digiData = reconData = 0 ;
    
    
    // external cel print
//    TTree * digiTree = cel.getTree(TestDigiComponent::name()) ;
//    digiTree->SetBranchAddress(TestDigiComponent::branchName(),&digiData) ;
//    TTree * reconTree = cel.getTree(TestReconComponent::name()) ;
//    reconTree->SetBranchAddress(TestDigiComponent::branchName(),&reconData) ;
//    std::cout << "Printing Events: " << std::endl ;
//    nEvents =  ;
//    for ( iEvent = 0 ; iEvent < lastEvt ; iEvent++ )
//     {
//      cel.deepRead(iEvent) ;
//      std::cout
//        <<"[main] Event "<<setw(2)<<std::right<<iEvent
//        <<" : "<<TestDigiComponent::name()<<" "<<*digiData
//        <<", "<<TestReconComponent::name()<<" "<<*reconData
//        <<std::endl ;
//     }
    delete digiData ;
    delete reconData ;
    digiData = reconData = 0 ;
    
    // print thanks to cel chains
    std::cout
      << "=================================\n"
      << "Printing events thanks to chains:\n"
      << std::flush ;
    delete digiData ;
    delete reconData ;
    digiData = reconData = 0 ;
    TObjArray * dataChains = new TObjArray ;
    TChain * celChain = cel.newChains(dataChains) ;
    TIter itrChain(dataChains) ;
    TChain * curChain ;
    while ( (curChain = (TChain*)itrChain.Next()) )
     {
      if (curChain->GetName()==TestDigiComponent::treeName())
       {
    	curChain->SetBranchAddress(TestDigiComponent::branchName(),&digiData) ;
        std::cout<<"Set address for Digi"<<std::endl ;
       }
      if (curChain->GetName()==TestReconComponent::treeName())
       {
    	curChain->SetBranchAddress(TestReconComponent::branchName(),&reconData) ;
        std::cout<<"Set address for Recon"<<std::endl ;
       }
     }    
    nEvents = celChain->GetEntries() ;
    for ( iEvent = 0 ; iEvent < nEvents ; iEvent++ )
     {
      celChain->LoadTree(iEvent) ;
      std::cout
        <<"[main] Event "<<setw(2)<<std::right<<iEvent
        <<" : "<<*digiData
        <<", "<<*reconData
        <<std::endl ;
     }
    delete digiData ;
    delete reconData ;
    digiData = reconData = 0 ;
    
    // print thanks to cel manager
    std::cout
      << "=====================================\n"
      << "Printing events thanks to CelManager:\n"
      << std::flush ;
    CelManager mgr(true) ;
    mgr.initRead(fileName) ;
    digiChain = mgr.getChainByType(TestDigiComponent::treeName()) ;
    digiChain->SetBranchAddress(TestDigiComponent::branchName(),&digiData) ;
    reconChain = mgr.getChainByType(TestReconComponent::treeName()) ;
    reconChain->SetBranchAddress(TestReconComponent::branchName(),&reconData) ;
    Long64_t iDigi, iRecon ;
    nEvents = mgr.getNumEvents() ;
    for ( iEvent = 0 ; iEvent < nEvents ; iEvent++ )
     {
      std::cout
        <<"[main] Event "<<setw(2)<<std::right<<iEvent
        <<" : "
        <<std::flush ;    	  
      iDigi = mgr.getEventIndexInTree(TestDigiComponent::treeName(),iEvent) ;
      digiChain->GetTree()->GetEntry(iDigi) ;
      iRecon = mgr.getEventIndexInTree(TestReconComponent::treeName(),iEvent) ;
      reconChain->GetTree()->GetEntry(iRecon) ;
      if ((digiData!=0)&&(reconData!=0))
       { std::cout<<*digiData<<", "<<*reconData<<std::endl ; }
      else
       { std::cout<<"Digi and/or Recon data not read"<<std::endl ; }
     }    
    digiChain = reconChain = 0 ;
    delete digiData ;
    delete reconData ;
    digiData = reconData = 0 ;
    
    
   }
  catch (...)
   {
    std::cout<<"AN UNKNOWN EXCEPTION HAS BEEN RAISED"<<std::endl ;
    sc = 1 ;
   }
         
  if (sc == 0)
   { std::cout << "TEST ROOT files reading succeeded!" << std::endl ; }
  else
   { std::cout << "FAILED recon reading" << std::endl ; }

  return(sc) ;
 }


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
    TestReader reader ;
    reader.add(baseName,new TestDigiComponent) ;
    reader.add(baseName,new TestReconComponent) ;
    reader.showByComponent() ;
    reader.showByEvent() ;
    
    CompositeEventList cel ;
    TString fileName = baseName ;
    fileName += ".cel.root" ;
    
    // internal cel print
    cel.openCelFile(fileName) ;
    cel.printInfo() ;
    cel.closeCelFile() ;
    
	// usual loop variables
    Long64_t nEvents, iEvent ;
    TestAbstractData * digiData =0 ;
    TestAbstractData * reconData =0 ;

    // external cel print
    delete digiData ;
    delete reconData ;
    digiData = reconData = 0 ;
    cel.openCelFile(fileName) ;
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
    cel.closeCelFile() ;
    
    // print thanks to cel chains
    std::cout << "Printing events thanks to chains: " << std::endl ;
    delete digiData ;
    delete reconData ;
    digiData = reconData = 0 ;
    cel.openCelFile(fileName) ;
    TObjArray * dataChains = new TObjArray ;
    TChain * celChain = cel.buildAllChains(dataChains) ;
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
    cel.closeCelFile() ;
    
    // print thanks to cel manager
    std::cout << "Printing events thanks to CelManager: " << std::endl ;
    delete digiData ;
    delete reconData ;
    digiData = reconData = 0 ;
    CelManager mgr(true) ;
    mgr.initRead(fileName) ;
    delete digiData ;
    digiData = 0 ;
    TChain * digiChain = mgr.getChainByType(TestDigiComponent::treeName()) ;
    digiChain->SetBranchAddress(TestDigiComponent::branchName(),&digiData) ;
    delete reconData ;
    reconData = 0 ;
    TChain * reconChain = mgr.getChainByType(TestReconComponent::treeName()) ;
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
    cel.closeCelFile() ;
    
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

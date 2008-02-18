
#include <rootUtil/TestData.h>
#include <rootUtil/CompositeEventList.h>
#include <rootUtil/CelManager.h>
#include <rootUtil/ComponentsInfoGlast.h>
#include <Riostream.h>

int main( int argc, char ** argv )
 {
  char * baseName = "test" ;
  int argi = 1 ;
  if (argc > 1)
   { baseName = argv[argi++] ; } 

  int sc = 0 ;
  try 
   {
    // description of components
    ComponentsInfoGlast infos ;
    const ComponentInfo * digiInfo = infos.getInfo("tdigi") ;
    const ComponentInfo * reconInfo = infos.getInfo("trecon") ;
    
    // direct reader access
    std::cout
      << "=================================\n"
      << "DIRECT READ\n"
      << std::flush ;
    TestReader reader ;
    reader.add(baseName,digiInfo) ;
    reader.add(baseName,reconInfo) ;
    reader.showByComponent() ;
    reader.showByEvent() ;
    

    // internal cel print
    std::cout
      << "=================================\n"
      << "CEL Internal meta info\n"
      << std::flush ;
    TString fileName = baseName ;
    fileName += ".cel.root" ;
    CompositeEventList cel(fileName) ;
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
    UInt_t digiIndex = cel.componentIndex(digiInfo->componentName) ;
    digiChain = cel.newChain(digiIndex) ;
    digiChain->SetBranchAddress(digiInfo->topBranchName,&digiData) ;
    UInt_t reconIndex = cel.componentIndex(reconInfo->componentName) ;
    reconChain = cel.newChain(reconIndex) ;
    reconChain->SetBranchAddress(reconInfo->topBranchName,&reconData) ;
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
    cel.setDataAddress(digiInfo->componentName,digiInfo->topBranchName,&digiData) ;
    cel.setDataAddress(reconInfo->componentName,reconInfo->topBranchName,&reconData) ;
    nEvents = cel.numEvents() ;
    for ( iEvent = 0 ; iEvent < nEvents ; iEvent++ )
     {
      std::cout
        <<"Event "<<setw(2)<<std::right<<iEvent
        <<" : "<<std::flush ;    	  
      cel.deepRead(iEvent) ;
      if ((digiData!=0)&&(reconData!=0))
       { std::cout<<*digiData<<", "<<*reconData<<std::endl ; }
      else
       { std::cout<<"Digi and/or Recon data not read"<<std::endl ; }
     }
    delete digiData ;
    delete reconData ;
    digiData = reconData = 0 ;
    
    
    // print thanks to cel chains
    std::cout
      << "=================================\n"
      << "FRIEND CHAIN\n"
      << std::flush ;
    delete digiData ;
    delete reconData ;
    digiData = reconData = 0 ;
    TObjArray * dataChains = new TObjArray ;
    TChain * celChain = cel.newChains(dataChains) ;
    TVirtualIndex * digiTreeIndex, * reconTreeIndex ;
    TIter itrChain(dataChains) ;
    TChain * curChain ;
    while ( (curChain = (TChain*)itrChain.Next()) )
     {
      if (curChain->GetName()==digiInfo->treeName)
       {
        digiChain = curChain ;
        digiTreeIndex = digiChain->GetTreeIndex() ;
    	curChain->SetBranchAddress(digiInfo->topBranchName,&digiData) ;
        std::cout<<"Set address for Digi"<<std::endl ;
       }
      if (curChain->GetName()==reconInfo->treeName)
       {
        reconChain = curChain ;
        reconTreeIndex = reconChain->GetTreeIndex() ;
    	curChain->SetBranchAddress(reconInfo->topBranchName,&reconData) ;
        std::cout<<"Set address for Recon"<<std::endl ;
       }
     }    
    nEvents = celChain->GetEntries() ;
    for ( iEvent = 0 ; iEvent < nEvents ; iEvent++ )
     {
      digiChain->SetTreeIndex(digiTreeIndex) ;
      reconChain->SetTreeIndex(reconTreeIndex) ;
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
      << "CelManager\n"
      << std::flush ;
    CelManager mgr(true) ;
    mgr.initRead(fileName) ;
    digiChain = mgr.getChainByType(digiInfo->treeName) ;
    digiChain->SetBranchAddress(digiInfo->topBranchName,&digiData) ;
    reconChain = mgr.getChainByType(reconInfo->treeName) ;
    reconChain->SetBranchAddress(reconInfo->topBranchName,&reconData) ;
    Long64_t iDigi, iRecon ;
    nEvents = mgr.getNumEvents() ;
    for ( iEvent = 0 ; iEvent < nEvents ; iEvent++ )
     {
      std::cout
        <<"[main] Event "<<setw(2)<<std::right<<iEvent
        <<" : "
        <<std::flush ;    	  
      iDigi = mgr.getEventIndexInTree(digiInfo->treeName,iEvent) ;
      digiChain->GetTree()->GetEntry(iDigi) ;
      iRecon = mgr.getEventIndexInTree(reconInfo->treeName,iEvent) ;
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

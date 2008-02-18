//
//========================================================
// File: ComponentsInfoGlast.cxx
// Main developer: D.Chamont.
//========================================================
//

#include <rootUtil/ComponentsInfoGlast.h>


ComponentsInfoGlast::ComponentsInfoGlast()
 {
  ComponentInfo * info ;
  
  // merit and auxiliary trees
  
  info = newInfo() ;
  info->componentName = "merit" ;
  info->treeName = "MeritTuple" ;
  info->runIdBranchName = "EvtRun" ;
  info->eventIdBranchName = "EvtEventId" ;
  
  info = newInfo() ;
  info->componentName = "pointing" ;
  info->treeName = "pointing_history" ;
  info->shortName = "pointing" ;
  
  info = newInfo() ;
  info->componentName = "jobinfo" ;
  info->treeName = "jobinfo" ;
  info->shortName = "jobinfo" ;
  
  // other tuple like trees
  
  info = newInfo() ;
  info->componentName = "merged" ;
  info->treeName = "MeritTuple" ;
  info->runIdBranchName = "EvtRun" ;
  info->eventIdBranchName = "EvtEventId" ;
  
  info = newInfo() ; // BTuple
  info->componentName = "beamtest" ; // BTuple
  info->treeName = "Output" ;
  info->runIdBranchName = "RunID" ;  
  info->eventIdBranchName = "EventID" ;  
  
  info = newInfo() ;
  info->componentName = "svac" ;
  info->treeName = "Output" ;
  info->runIdBranchName = "RunID" ;  
  info->eventIdBranchName = "EventID" ;  
  
  info = newInfo() ;
  info->componentName = "cal" ;
  info->treeName = "CalTuple" ;
  info->runIdBranchName = "RunID" ;  
  info->eventIdBranchName = "EventID" ;  
  
  // deep trees 
 
  info = newInfo() ;
  info->componentName = "digi" ;
  info->treeName = "Digi" ;
  info->topBranchName = "DigiEvent" ;
  info->topBranchType = "DigiEvent" ;
  info->libName = "libdigiRootData.so" ;
  info->runIdBranchName = "m_runId" ;
  info->eventIdBranchName = "m_eventId" ;
  
  info = newInfo() ;
  info->componentName = "recon" ;
  info->treeName = "Recon" ;
  info->topBranchName = "ReconEvent" ;
  info->topBranchType = "ReconEvent" ;
  info->libName = "libreconRootData.so" ;
  info->runIdBranchName = "m_runId" ;
  info->eventIdBranchName = "m_eventId" ;
  
  info = newInfo() ;
  info->componentName = "mc" ;
  info->treeName = "Mc" ;
  info->topBranchName = "McEvent" ;
  info->topBranchType = "McEvent" ;
  info->libName = "libmcRootData.so" ;
  info->runIdBranchName = "m_runId" ;
  info->eventIdBranchName = "m_eventId" ;
  
  info = newInfo() ;
  info->componentName = "gcr" ;
  info->treeName = "GcrSelect" ;
  info->topBranchName = "GcrSelectEvent" ;
  info->topBranchType = "GcrSelectEvent" ;
  info->libName = "libgcrSelectRootData.so" ;
  info->runIdBranchName = "m_runId" ;
  info->eventIdBranchName = "m_eventId" ;
  
  // for rootUtil tests
  
  // Below, one cannot use the typedefs for the topBranchType
  // values, because it is reused for TTree:Branch(), which
  // does not know about typedefs.

  info = newInfo() ;
  info->componentName = "tdigi" ;
  info->treeName = "TestDigiTree" ;
  info->topBranchName = "TestDigiBranch" ;
  info->topBranchType = "TestData<TestDigiLabel>" ;
  
  info = newInfo() ;
  info->componentName = "trecon" ;
  info->treeName = "TestReconTree" ;
  info->topBranchName = "TestReconBranch" ;
  info->topBranchType = "TestData<TestReconLabel>" ;

 }
 

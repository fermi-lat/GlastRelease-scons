#define McSelector_cxx
// The class definition in McSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector().
//
// This class is derived from the ROOT class TSelector.
// The following members functions are called by the TTree::Process() functions:
//    Begin():       called everytime a loop on the tree starts,
//                   a convenient place to create your histograms.
//    Notify():      this function is called at the first entry of a new Tree
//                   in a chain.
//    Process():     called for each event. In this function you decide what 
//                   to read and you fill your histograms.
//    Terminate():   called at the end of a loop on the tree,
//                   a convenient place to draw/fit your histograms.
//
//   To use this file, try the following session on your Tree T
//
// Root > T->Process("McSelector.C")
// Root > T->Process("McSelector.C","some options")
// Root > T->Process("McSelector.C+")
//
#include "McSelector.h"
#include "TH2.h"
#include "TStyle.h"

ClassImp(McSelector)

void McSelector::Begin(TTree *tree)
{
   // Function called before starting the event loop.
   // When running with PROOF Begin() is only called in the client.
   // Initialize the tree branches.


   TString option = GetOption();

   // Create the /root/glast folder, if necessary
   glastFolder = (TFolder*)gROOT->FindObjectAny("/glast");
   if (!glastFolder)
       glastFolder = gROOT->GetRootFolder()->AddFolder("glast","glast top level folders");
   mcFolder = (TFolder*)gROOT->FindObjectAny("/glast/mc");
   if (!mcFolder)
       mcFolder = glastFolder->AddFolder("mc","monte carlo folder");
   histoFolder = (TFolder*)gROOT->FindObjectAny("/glast/histograms");
   if (!histoFolder) 
       histoFolder = glastFolder->AddFolder("histograms", "histogram folder");
   taskFolder = (TFolder*)gROOT->FindObjectAny("/glast/tasks");
   if (!taskFolder) 
       std::cout << "No task folder - no tasks to process!" << std::endl;
   else {
       m_mcTask = (GlastMcTask*)taskFolder->FindObjectAny("GlastMcTask");
       if (!m_mcTask) {
           std::cout << "NO TASKS FOUND" << std::endl;
       } else {
           m_mcTask->CreateHistograms();
       }
   }

  std::cout <<"done with begin of selector" << std::endl;
}

void McSelector::SlaveBegin(TTree *tree)
{
   // Function called before starting the event loop.
   // When running with PROOF SlaveBegin() is called in each slave
   // Initialize the tree branches.

   Init(tree);

   TString option = GetOption();

}

Bool_t McSelector::Process(Int_t entry)
{
   // Processing function. This function is called
   // to process an event. It is the user's responsability to read
   // the corresponding entry in memory (may be just a partial read).
   // Once the entry is in memory one can apply a selection and if the
   // event is selected histograms can be filled.

   if (!m_mcTask) {
     std::cout << "NO TASKS" << std::endl;
     return kTRUE;
   }

   std::cout <<"requesting data entry = " << entry << std::endl;
   RequestedMcData *rmd = (RequestedMcData*)taskFolder->FindObjectAny("RequestedMcData");
   if(rmd && rmd->McPartData()) {
     fChain->GetEntry(entry);
     mcFolder->Add(m_mcEvent);
     //b_McEvent_m_mcParticleCol->GetEntry(entry);
     //mcFolder->Add(m_eCol);
   } else {
     std::cout << "McParticle Col not requested" << std::endl;
   }

   std::cout << "calling exe" << std::endl;
   m_mcTask->ExecuteTask();

   return kTRUE;

}

void McSelector::SlaveTerminate()
{
   // Function called at the end of the event loop in each PROOF slave.


}

void McSelector::Terminate()
{
   // Function called at the end of the event loop.
   // When running with PROOF Terminate() is only called in the client.


}

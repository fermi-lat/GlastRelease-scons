#define DigiSelector_cxx
// The class definition in DigiSelector.h has been generated automatically
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
// Root > T->Process("DigiSelector.C")
// Root > T->Process("DigiSelector.C","some options")
// Root > T->Process("DigiSelector.C+")
//
#include "DigiSelector.h"
#include "TH2.h"
#include "TStyle.h"

ClassImp(DigiSelector)

void DigiSelector::Begin(TTree *tree)
{
   // Function called before starting the event loop.
   // When running with PROOF Begin() is only called in the client.
   // Initialize the tree branches.

   TString option = GetOption();
   // Create the /root/glast folder, if necessary
   glastFolder = (TFolder*)gROOT->FindObjectAny("/glast");
   if (!glastFolder)
       glastFolder = gROOT->GetRootFolder()->AddFolder("glast","glast top level folders");
   digiFolder = (TFolder*)gROOT->FindObjectAny("/glast/digi");
   if (!digiFolder)
       digiFolder = glastFolder->AddFolder("digi","digitization folder");
   histoFolder = (TFolder*)gROOT->FindObjectAny("/glast/histograms");
   if (!histoFolder) 
       histoFolder = glastFolder->AddFolder("histograms", "histogram folder");
   taskFolder = (TFolder*)gROOT->FindObjectAny("/glast/tasks");
   if (!taskFolder) 
       std::cout << "No task folder - no tasks to process!" << std::endl;
   else {
       m_digiTask = (GlastDigiTask*)taskFolder->FindObjectAny("DigiTask");
       if (!m_digiTask) {
           std::cout << "NO TASKS FOUND" << std::endl;
       } else {
           m_digiTask->CreateHistograms();
       }
   }
}

void DigiSelector::SlaveBegin(TTree *tree)
{
   // Function called before starting the event loop.
   // When running with PROOF SlaveBegin() is called in each slave
   // Initialize the tree branches.

   Init(tree);

   TString option = GetOption();

}

Bool_t DigiSelector::Process(Int_t entry)
{
   // Processing function. This function is called
   // to process an event. It is the user's responsability to read
   // the corresponding entry in memory (may be just a partial read).
   // Once the entry is in memory one can apply a selection and if the
   // event is selected histograms can be filled.

   // WARNING when a selector is used with a TChain, you must use
   //  the pointer to the current Tree to call GetEntry(entry).
   //  entry is always the local entry number in the current tree.
   //  Assuming that fChain is the pointer to the TChain being processed,
   //  use fChain->GetTree()->GetEntry(entry);
    //fChain->GetTree()->GetEntry(entry);
    if (!m_digiTask) {
        std::cout <<"NOTASKS" << std::endl;
        return kTRUE;
    }

   RequestedDigiData *rdd = (RequestedDigiData*)taskFolder->FindObjectAny("RequestedDigiData");
   if (rdd && rdd->AcdData()) {
        b_DigiEvent_acdDigiCol->GetEntry(entry);
        digiFolder->Add(m_acdDigiCol);
    } else {
        std::cout << "ACD not requested" << std::endl;
    }
    m_digiTask->ExecuteTask();

   return kTRUE;
}

void DigiSelector::SlaveTerminate()
{
   // Function called at the end of the event loop in each PROOF slave.


}

void DigiSelector::Terminate()
{
   // Function called at the end of the event loop.
   // When running with PROOF Terminate() is only called in the client.


}

//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Thu May 13 22:36:25 2004 by ROOT version3.10/02)
//   from TTree Recon/GLAST Reconstruction Data
//   found on file: recon_1500.root
//////////////////////////////////////////////////////////


#ifndef ReconSelector_h
#define ReconSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TFolder.h>
#if !defined(__CINT__)
#include <iostream>
#else  // for interactive use
#include "Riostream.h"
#endif
#include "reconRootData/ReconEvent.h"
class ReconEvent;
#include "GlastReconTask.h"

class ReconSelector : public TSelector {
   public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TFolder        *glastFolder, *reconFolder, *histoFolder, *taskFolder;
   GlastReconTask *m_reconTask;

//Declaration of leaves types
   ReconEvent      *m_reconEvent;
   UInt_t          m_eventId;
   UInt_t          m_runId;
   AcdRecon        *m_acd;
   CalRecon        *m_cal;
   TkrRecon        *m_tkr;

//List of branches
   TBranch        *b_ReconEvent_m_eventId;   //!
   TBranch        *b_ReconEvent_m_runId;   //!
   TBranch        *b_ReconEvent_m_acd; //!
   TBranch        *b_ReconEvent_m_cal; //!
   TBranch        *b_ReconEvent_m_tkr;  //!
 
   ReconSelector(TTree *tree=0) { }
   ~ReconSelector() { }
   Int_t   Version() const {return 1;}
   void    Begin(TTree *tree);
   void    SlaveBegin(TTree *tree);
   void    Init(TTree *tree);
   Bool_t  Notify();
   Bool_t  Process(Int_t entry);
   void    SetOption(const char *option) { fOption = option; }
   void    SetObject(TObject *obj) { fObject = obj; }
   void    SetInputList(TList *input) {fInput = input;}
   TList  *GetOutputList() const { return fOutput; }
   void    SlaveTerminate();
   void    Terminate();
   ClassDef(ReconSelector,0)
};

#endif

#ifdef ReconSelector_cxx
void ReconSelector::Init(TTree *tree)
{
//   Set branch addresses
   if (tree == 0) return;
   fChain    = tree;
   // Comment this out for now
   //fChain->SetMakeClass(1);

   m_reconEvent =0;
   fChain->SetBranchAddress("ReconEvent", &m_reconEvent);
/*
   fChain->SetBranchAddress("m_eventId",&m_eventId);
   fChain->SetBranchAddress("m_runId",&m_runId);
   m_acd =  new AcdRecon();
   m_cal = new CalRecon();
   m_tkr = new TkrRecon();

   fChain->SetBranchAddress("m_acd", &m_acd);
   fChain->SetBranchAddress("m_cal", &m_cal);
   fChain->SetBranchAddress("m_tkr", &m_tkr);
*/
}

Bool_t ReconSelector::Notify()
{
   // Called when loading a new file.
   // Get branch pointers.
/*   b_ReconEvent_m_eventId = fChain->GetBranch("m_eventId");
   b_ReconEvent_m_runId = fChain->GetBranch("m_runId");
   b_ReconEvent_m_acd = fChain->GetBranch("m_acd");
   b_ReconEvent_m_cal = fChain->GetBranch("m_cal");
   b_ReconEvent_m_tkr = fChain->GetBranch("m_tkr");
*/
   return kTRUE;
}

#endif // #ifdef ReconSelector_cxx


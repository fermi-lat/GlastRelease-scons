//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Thu May 13 22:35:27 2004 by ROOT version3.10/02)
//   from TTree Mc/GLAST Monte Carlo Data
//   found on file: mc_1500.root
//////////////////////////////////////////////////////////


#ifndef McSelector_h
#define McSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TFolder.h>
#include "mcRootData/McEvent.h"
#if !defined(__CINT__)
#include <iostream>
#else  // for interactive use
#include "Riostream.h"
#endif
#include "GlastMcTask.h"
class McEvent;

class McSelector : public TSelector {
   public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TFolder        *glastFolder, *mcFolder, *histoFolder, *taskFolder;
//Declaration of leaves types
   GlastMcTask     *m_mcTask;
   McEvent         *m_mcEvent;
   UInt_t          m_eventId;
   UInt_t          m_runId;
   Int_t           m_sourceId;
   UInt_t          m_sequence;
   Double_t        m_timeStamp;
   TObjArray       *m_mcParticleCol;
   TObjArray       *m_mcIntHitCol;
   TObjArray       *m_mcPosHitCol;

//List of branches
   TBranch        *b_McEvent_m_eventId;   //!
   TBranch        *b_McEvent_m_runId;   //!
   TBranch        *b_McEvent_m_sourceId;   //!
   TBranch        *b_McEvent_m_sequence;   //!
   TBranch        *b_McEvent_m_timeStamp;   //!
   TBranch        *b_McEvent_m_mcParticleCol; //!
   TBranch        *b_McEvent_m_mcIntHitCol; //!
   TBranch        *b_McEvent_m_mcPosHitCol; //!

   McSelector(TTree *tree=0) { }
   ~McSelector() { }
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
   ClassDef(McSelector,0)
};

#endif

#ifdef McSelector_cxx
void McSelector::Init(TTree *tree)
{
//   Set branch addresses
   if (tree == 0) return;
   fChain    = tree;
   // Comment this out for now
   //fChain->SetMakeClass(1);

   m_mcEvent = 0;
   fChain->SetBranchAddress("McEvent", &m_mcEvent);

/* Seem to have only one branch right now.. need to load whole McEvent
   m_mcParticleCol = new TObjArray();
   m_mcIntHitCol = new TObjArray();
   m_mcPosHitCol = new TObjArray();

   fChain->SetBranchAddress("m_eventId",&m_eventId);
   fChain->SetBranchAddress("m_runId",&m_runId);
   fChain->SetBranchAddress("m_sourceId",&m_sourceId);
   fChain->SetBranchAddress("m_sequence",&m_sequence);
   fChain->SetBranchAddress("m_timeStamp",&m_timeStamp);
 
   fChain->SetBranchAddress("m_particleCol",&m_mcParticleCol);
   fChain->SetBranchAddress("m_integratingHitCol", &m_mcIntHitCol);
   fChain->SetBranchAddress("m_positionHitCol", &m_mcPosHitCol);
 */
}

Bool_t McSelector::Notify()
{
   // Called when loading a new file.
   // Get branch pointers.
/*
   b_McEvent_m_eventId = fChain->GetBranch("m_eventId");
   b_McEvent_m_runId = fChain->GetBranch("m_runId");
   b_McEvent_m_sourceId = fChain->GetBranch("m_sourceId");
   b_McEvent_m_sequence = fChain->GetBranch("m_sequence");
   b_McEvent_m_timeStamp = fChain->GetBranch("m_timeStamp");

   b_McEvent_m_mcParticleCol = fChain->GetBranch("m_particleCol");
   b_McEvent_m_mcIntHitCol = fChain->GetBranch("m_integratingHitCol");
   b_McEvent_m_mcPosHitCol = fChain->GetBranch("m_positionHitCol");
*/
   return kTRUE;
}

#endif // #ifdef McSelector_cxx


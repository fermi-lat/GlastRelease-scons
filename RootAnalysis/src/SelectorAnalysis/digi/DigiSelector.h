//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Fri Apr 30 23:17:18 2004 by ROOT version3.10/02)
//   from TTree Digi/Digi
//   found on file: gsfc_acd_test_digi.root
//////////////////////////////////////////////////////////


#ifndef DigiSelector_h
#define DigiSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#if !defined(__CINT__)
#include <iostream>
#else  // for interactive use
#include "Riostream.h"
#endif
#include "digiRootData/DigiEvent.h"
class DigiEvent;
#include "GlastDigiTask.h"

   const Int_t kMaxm_acdDigiCol = 20;
   const Int_t kMaxm_calDigiCloneCol = 1;

class DigiSelector : public TSelector {
   public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TFolder        *glastFolder, *digiFolder, *histoFolder, *taskFolder;
   GlastDigiTask        *m_digiTask;
//Declaration of leaves types
   DigiEvent       *m_digiEvent;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Double_t        m_timeStamp;
   UInt_t          m_runId;
   UInt_t          m_eventId;
   UChar_t         m_fromMc;
   UInt_t          m_levelOneTrigger_fUniqueID;
   UInt_t          m_levelOneTrigger_fBits;
   UInt_t          m_levelOneTrigger_m_trigger;

   TClonesArray   *m_acdDigiCol;
   Int_t           m_numAcdDigis;
   TClonesArray   *m_calDigiCol;
   Int_t           m_numCalDigis;
   TObjArray       *m_tkrDigiCol;

//List of branches
   TBranch        *b_DigiEvent_fUniqueID;   //!
   TBranch        *b_DigiEvent_fBits;   //!
   TBranch        *b_DigiEvent_timeStamp;   //!
   TBranch        *b_DigiEvent_runId;   //!
   TBranch        *b_DigiEvent_eventId;   //!
   TBranch        *b_DigiEvent_fromMc;   //!
   TBranch        *b_DigiEvent_levelOneTrigger;//!
   TBranch        *b_DigiEvent_acdDigiCol; //!
   TBranch        *b_DigiEvent_numAcdDigis;   //!
   TBranch        *b_DigiEvent_calDigiCol; //!
   TBranch        *b_DigiEvent_numCalDigis;   //!
   TBranch        *b_DigiEvent_tkrDigiCol; //!

   DigiSelector(TTree *tree=0) { }
   ~DigiSelector() { }
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
public:
   ClassDef(DigiSelector,0)
};

#endif

#ifdef DigiSelector_cxx
void DigiSelector::Init(TTree *tree)
{

//   Set branch addresses
   if (tree == 0) return;
   fChain    = tree;
   // Use this if we are accessing individual branches rather than whole digiEvent
   // SetMakeClass doesn't work if we already loaded the library - it seems
   fChain->SetMakeClass(1);

   m_acdDigiCol = new TClonesArray("AcdDigi", 1);
   fChain->SetBranchAddress("m_acdDigiCol", &m_acdDigiCol);
   fChain->SetBranchAddress("m_numAcdDigis",&m_numAcdDigis);
   m_calDigiCol = new TClonesArray("CalDigi", 1);
   fChain->SetBranchAddress("m_calDigiCol", &m_calDigiCol);

}

Bool_t DigiSelector::Notify()
{
   // Called when loading a new file.
   // Get branch pointers.


   b_DigiEvent_acdDigiCol = fChain->GetBranch("m_acdDigiCol");
   b_DigiEvent_numAcdDigis = fChain->GetBranch("m_numAcdDigis");
   b_DigiEvent_calDigiCol = fChain->GetBranch("m_calDigiCol");

   return kTRUE;
}


#endif // #ifdef DigiSelector_cxx


#include "RootTreeAnalysis.h"
#include "TObjArray.h"

UInt_t digiEventId, reconEventId, mcEventId;
UInt_t digiRunNum, reconRunNum, mcRunNum;


#include "McFragment.cxx"


/* Process the Digi Data
   Called by Go()
*/

#include "DigiCalFragment.cxx"

#include "TrgFragment.cxx"

#include "DigiAcdFragment.cxx"
#include "ReconAcdFragment.cxx"


#include "ReconCalFragment.cxx"
#include "ReconTkrFragment.cxx"

void RootTreeAnalysis::DigiTkrHistDefine() {};
void RootTreeAnalysis::DigiTkr() {};

/* Event Loop
All Analysis goes here
*/
void RootTreeAnalysis::Go(Int_t numEvents)
{    
    //  To read only selected branches - saves processing time
    //  Comment out any branches you are not interested in.

    // mc branches:
    if (mcTree) {
        mcTree->SetBranchStatus("*", 0);    // disable all branches
        // Activate desired branches...
        mcTree->SetBranchStatus("m_eventId", 1);
        mcTree->SetBranchStatus("m_particleCol", 1);
        mcTree->SetBranchStatus("m_runId", 1);        
        mcTree->SetBranchStatus("m_integratingHitCol", 1);        
        mcTree->SetBranchStatus("m_positionHitCol", 1);        
    }
    
    
    // determine how many events to process
    Int_t nentries = GetEntries();
    std::cout << "\nNum Events in File is: " << nentries << std::endl;
    Int_t curI;
    Int_t nMax = TMath::Min(numEvents+m_StartEvent,nentries);
    if (m_StartEvent == nentries) {
        std::cout << " all events in file read" << std::endl;
        return;
    }
    if (nentries <= 0) return;
    
    // Keep track of how many bytes we have read in from the data files
    Int_t nbytes = 0, nb = 0;

    // BEGINNING OF EVENT LOOP
    for (Int_t ievent=m_StartEvent; ievent<nMax; ievent++, curI=ievent) {
        
      // Too verbose
      // mc->Dump();

        if (mc) {
          mc->Clear();
        }
        
        if (evt) {
          evt->Clear();
        }
        
        if (rec) {
          rec->Clear();
        }
        
        digiEventId = 0; reconEventId = 0; mcEventId = 0;
        digiRunNum = 0; reconRunNum = 0; mcRunNum = 0;
        
        nb = GetEvent(ievent);
        nbytes += nb;
        
        
        // Monte Carlo ONLY analysis
        if (mc) {  // if we have mc data process it
            mcEventId = mc->getEventId();
            mcRunNum = mc->getRunId();
            McData();
        } 
        
        // Digi ONLY analysis
        if (evt) {
            digiEventId = evt->getEventId(); 
            digiRunNum = evt->getRunId();            
	    // DigiAcd();
	    // DigiTkr();
            DigiCal();
	    Trigger();
        }
        
        // Recon ONLY analysis
        if (rec) {  // If we have mc data process it
            reconEventId = rec->getEventId();
            reconRunNum = rec->getRunId();
	    // ReconAcd();
	    // ReconTkr();
            ReconCal();
        } 
        
    }  // end analysis code in event loop
    
    m_StartEvent = curI;

    return;
}

/*! Digi Analysis Histogram Defintions */

void RootTreeAnalysis::DigiHistDefine() {

  histFile->cd();

  DigiCalHistDefine();
  // DigiTkrHistDefine();
  // DigiAcdHistDefine();
  TriggerHistDefine();

}

/*! Recon Analysis Histogram Defintions */

void RootTreeAnalysis::ReconHistDefine() {

  histFile->cd();

  ReconCalHistDefine();
  //  ReconTkrHistDefine();
  //  ReconAcdHistDefine();

};    


/* Setup ALL histograms */
void RootTreeAnalysis::HistDefine() {
    
    gStyle->SetOptStat(111111);
    
    histFile = new TFile(m_histFileName,"RECREATE");
    
    McHistDefine();
    DigiHistDefine();
    ReconHistDefine();
}


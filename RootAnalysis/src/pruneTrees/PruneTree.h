/** @class RootTreeAnalysis
* @brief   This class can manipulate a Monte Carlo, Digi, and a Recon Root file
* at the same time.  It is meant as a full example of using ROOT to manipulate
* GLAST ROOT data.  
*
* This class is intended to provide useful manipulation
* of the Root Event loop. Users put their analysis code
* into the Go function (in RootTreeAnalysis.cxx). They should not
* need to look at RootTreeAnalysis.h (except to see the interface)
* 
* allows for:
* init and re-init use of Root file(s)
* clear all histograms
* 'go n events' allowing to continue on in the file
* or 'rewind'
* 
*  Example of use:
*
*  gROOT->LoadMacro("RootTreeAnalysis.cxx");     // 'compile' class
*  RootTreeAnalysis m("myDigiFile.root", "myReconFile.root"); // create RootTreeAnalysis object
*  m.Go(500);      // loop over 500 events. Go contains your analysis code
*  ... look at histograms ...
*  m.Go()          // look at remainder of file
*  ... look at histograms ...
*  m.HistClear();      // clear histograms
*  m.Init("AnotherRootFile.root");
*  m.Go(50);
*  ... and so on ...
*
* After editing your Go function, you need to issue a gROOT->Reset() and
* repeat the above sequence starting from the .L RootTreeAnalysis.cxx.
*
* If you only have a digi or only a recon root file... setup RootTreeAnalysis like this:
* if you only have a recon root file:
* RootTreeAnalysis *m = new RootTreeAnalysis("", "myReconFile.root")
* if you only have a digi root file:
* RootTreeAnalysis *m = new RootTreeAnalysis("myDigiFile.root", "")
*
* Version 0.1 17-Mar-1999 Richard Creation
* Version 1.0 Spring, 2000 Revised for use with GLAST 1999 TestBeam
* Version 1.5 25-Mar-2001  Revised for use with GLAST 2001 Balloon
* Version 2.0 14-Aug-2001 Final version for GLAST 2001 Balloon
*/


#ifndef PruneTree_h
#define PruneTree_h 1

#if !defined(__CINT__)
// Need these includes if we wish to compile this code
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h"
#include "TStyle.h"
#include "TCollection.h"  // Declares TIter
#include "digiRootData/DigiEvent.h"
#include "reconRootData/ReconEvent.h"
#include "mcRootData/McEvent.h"
#include <iostream>
#else  // for interactive use
#include "iostream.h"
class DigiEvent;
class ReconEvent;
class McEvent;
#endif

class RootTreeAnalysis {
public :
    /// Histogram file
    TH1F        phaArr_A[500];
    TH1F        phaArr_B[500];
    TFile       *histFile;
    /// Input digitization file
    TFile       *digiFile;   
    /// Input reconstruction file
    TFile       *reconFile;  
    /// Input monte carlo file
    TFile       *mcFile;     
    /// Input ntuple file
    TFile       *anaTupFile;     
    /// pointer to the digi tree
    TTree       *digiTree;    
    /// pointer to the reconstruction tree
    TTree       *reconTree;
    /// pointer to the ntuple tree
    TTree       *anaTupTree;
    /// pointer to the monte carlo tree
    TTree       *mcTree;      
    /// Copy Input digitization file
    TFile       *copyDigiFile;   
    /// Copy Input reconstruction file
    TFile       *copyReconFile;  
    /// Copy Input monte carlo file
    TFile       *copyMcFile;     
    /// Copy Input ntuple file
    TFile       *copyAnaTupFile;     
    /// pointer to the digi tree Copy
    TTree       *copyDigiTree;    
    /// pointer to the reconstruction tree Copy
    TTree       *copyReconTree;
    /// pointer to the ntuple tree Copy
    TTree       *copyAnaTupTree;
    /// pointer to the monte carlo tree Copy
    TTree       *copyMcTree;      
    /// Optional TChain input
    TChain      *m_digiChain, *m_recChain, *m_mcChain, *m_anaTupChain;
    /// Optional TChain input Copy
    TChain      *m_copyDigiChain, *m_copyRecChain, *m_copyMcChain, *m_copyAnaTupChain;
    /// pointer to a DigiEvent
    DigiEvent   *evt;
    /// pointer to a ReconEvent
    ReconEvent  *rec;
    /// Pointer to a McEvent
    McEvent     *mc;
    /// array of doubles for ntuple branches
    double     anaTup[1000];
    /// name of the output histogram ROOT file
    char        *m_histFileName; 
    /// Arrays that contain pointers to the TFile, TTree, and TChains
    TObjArray   *fileArr, *treeArr, *chainArr;
    /// Arrays that contain pointers to the copied TFile, TTree, and TChains
    TObjArray   *fileArrCopy, *treeArrCopy, *chainArrCopy;
    
    int pruned;

	/// Default ctor, requires that that user calls RootTreeAnalysis::Init
	/// to setup access to specific ROOT files.
    RootTreeAnalysis(); 
    
	/// Standard ctor, where user provides the names of the input root files
	/// and optionally the name of the output ROOT histogram file
    RootTreeAnalysis( 
        const char* digiFileName, 
        const char* reconFileName="", 
        const char* mcFileName="", 
        const char* anaTupFileName="", 
        char *histFileName="Histograms.root"); 

	/// Special ctor which accepts TChains for input files
    RootTreeAnalysis( 
        TChain *digiChain, 
        TChain *recChain = 0, 
        TChain *mcChain = 0, 
        TChain *anaTupChain = 0, 
        char *histFileName="Histograms.root");

    ~RootTreeAnalysis();  

    /// start next Go with this event
    void StartWithEvent(Long64_t event) { m_StartEvent = event; };  
    /// reset for next Go to start at beginning of file 
    void Rewind() { m_StartEvent = 0; }; 
    /// allow the user to specify their own file name for the output ROOT file
    void SetHistFileName(char *histFileName) { m_histFileName = histFileName; }
    /// re-init with these input Root files
    void Init(  const char* digiFileName="", 
		const char* reconFileName="", 
		const char* mcFileName="", 
		const char* anaTupFileName="");
    /// define user histograms, ntuples and other output objects that will be saved to output
    void HistDefine();   
    /// make list of user histograms and all objects created for output
    void MakeHistList(); 
    /// write the existing histograms and ntuples out to file
    void WriteHist() { if (histFile) histFile->Write(); }; 
    /// Reset() all user histograms
    void HistClear(); 
    /// Retrieve a pointer to an object stored in our output ROOT file
    TObject* GetObjectPtr(const char *tag) { return (m_histList->FindObject(tag)); };
    /// process events
    void Go(Long64_t numEvents=100000000); 
    /// returns number of events in all open files
    Long64_t GetEntries() const;
    /// retrieve a pointer to event number.
    UInt_t GetEvent(Long64_t ievt);
    UInt_t CopyTrees();
    UInt_t CopyEntry(Long64_t ievt);
    
private:
    /// starting event number
    Long64_t m_StartEvent;
    /// list of user histograms
    THashList *m_histList;
        
    /// reset all member variables
    void Clear(); 

    /// Setup the Monte Calro output histograms
    void McHistDefine();
	/// Setup the Digitization output histograms
    void DigiHistDefine();
	/// Setup the Reconstruction output histograms
    void ReconHistDefine();
    
    /// event processing for the monte carlo data
    void McData();

    /// event processing for the digi TKR data
    void DigiTkr();
	/// event processing for digi CAL data
    void DigiCal();
	/// event processing for digi ACD data
    void DigiAcd();
    
    /// event processing for the recon TKR data
    void ReconTkr();
	/// event processing for the recon CAL data
    void ReconCal();
	/// event processing for the recon ACD data
    void ReconAcd();
    
};


inline RootTreeAnalysis::RootTreeAnalysis() 
{
    Clear();
}

inline RootTreeAnalysis::RootTreeAnalysis(const char* digiFileName, 
                                   const char* reconFileName, 
                                   const char* mcFileName, 
                                   const char* anaTupFileName, 
                                   char* histFileName)
{
	// Purpose and Method:  Standard constructor where the user provides the 
	//  names of input ROOT files and optionally the name of the output ROOT
	//  histogram file.
    printf(" opening files:\n\tdigi:\t%s\n\trecon:\t%s\n\tmc:\t%s\n\tanaTup:%s\n",
		digiFileName, reconFileName, mcFileName, anaTupFileName);
    
    Clear();
    
    SetHistFileName(histFileName);
    HistDefine();
    MakeHistList();
    
    Init(digiFileName, reconFileName, mcFileName, anaTupFileName);
}

inline RootTreeAnalysis::RootTreeAnalysis(TChain *digiChain, 
                                   TChain *recChain, 
                                   TChain *mcChain, 
				   TChain *anaTupChain,
                                   char *histFileName) {


    Clear();
    
    SetHistFileName(histFileName);
    HistDefine();
    MakeHistList();
    
    if (chainArr) delete chainArr;
    chainArr = new TObjArray();
    
    if (mcChain != 0) {
        m_mcChain = mcChain;
        mc = 0;
        m_mcChain->SetBranchAddress("McEvent",&mc);
        chainArr->Add(m_mcChain);
    }

    if (digiChain != 0) {
        evt = 0;
        m_digiChain = digiChain;
        m_digiChain->SetBranchAddress("DigiEvent",&evt);
        chainArr->Add(m_digiChain);
    }
    
    if (recChain != 0) {
        m_recChain = recChain;
        rec = 0;
        m_recChain->SetBranchAddress("ReconEvent",&rec);
        chainArr->Add(m_recChain);
    }

    if (anaTupChain != 0) {
        m_anaTupChain = anaTupChain;
	TObjArray* bList = m_anaTupChain->GetListOfBranches();
	for (Long64_t ib=0; ib<bList->GetEntries(); ib++) {
	  const char* name = bList->At(ib)->GetName();
	  m_anaTupChain->SetBranchAddress(name,&anaTup[ib]);
	}
        chainArr->Add(m_anaTupChain);
    }
    m_StartEvent = 0;
    pruned = 0;
}


inline RootTreeAnalysis::~RootTreeAnalysis() {
    histFile->Close();
    
    //if (m_histList) delete m_histList;
    
    if (histFile) delete histFile;
    
    if (digiFile) delete digiFile;
    if (reconFile) delete reconFile;
    if (mcFile) delete mcFile;
    
    if (evt) { 
		evt->Clear(); 
		delete evt;
	}
    if (rec) {
		rec->Clear();
		delete rec;
	}
    if (mc) {
		mc->Clear();
		delete mc;
	}
    
    digiTree = 0;
    reconTree = 0;
    mcTree = 0;
    anaTupTree = 0;
    
    if (fileArr) delete fileArr;
    if (treeArr) delete treeArr;
    if (chainArr) delete chainArr;

    if (treeArrCopy) {
        for (Long64_t i = 0; i < treeArrCopy->GetEntries(); i++) {
            ((TTree*)treeArrCopy->At(i))->AutoSave();
            ((TTree*)treeArrCopy->At(i))->Print();
        }
    }
    if (chainArrCopy) {
        for (Long64_t i = 0; i < chainArrCopy->GetEntries(); i++) {
            ((TChain*)chainArrCopy->At(i))->AutoSave();
            ((TChain*)chainArrCopy->At(i))->Print();
        }
    }

    if (treeArrCopy) delete treeArrCopy;
    if (fileArrCopy) delete fileArrCopy;
    if (chainArrCopy) delete chainArrCopy;

    std::cout << "Events surviving pruning " << pruned << std::endl;

	Clear();
}


inline void RootTreeAnalysis::Init(const char* digiFileName, const char* reconFileName, const char* mcFileName, const char* anaTupFileName)
{
    // Purpose and Method:  Re-initialize file, tree, event pointers, using the 
	//   input ROOT files.  Histograms are *not* cleared.
    
    if (fileArr) delete fileArr;
    fileArr = new TObjArray();
    
    if (treeArr) delete treeArr;
    treeArr = new TObjArray();
         
    if (mcFile) {
        delete mc; 
        mc = 0;
        mcTree = 0;
        delete mcFile; 
        mcFile = 0;
    }
    
    if (mcFileName != "") {
        mcFile = new TFile(mcFileName);
        if (mcFile->IsOpen() == kTRUE) {
            mcTree = (TTree*)gDirectory->Get("Mc");
            mc = 0;
            mcTree->SetBranchAddress("McEvent",&mc);
            fileArr->Add(mcFile);
            treeArr->Add(mcTree);
        } else {
            mcFile = 0;
            std::cout << "mc data file could not be opened!!" << std::endl;
        }
    }

    if (digiFile) {
        delete evt; 
        evt = 0;
        digiTree = 0;
        delete digiFile; 
        digiFile = 0;
    }
    
    if (digiFileName != "") {
        digiFile = new TFile(digiFileName);
        if (digiFile->IsOpen() == kTRUE) {
            digiTree = (TTree*)gDirectory->Get("Digi");
            evt = 0;
            digiTree->SetBranchAddress("DigiEvent",&evt);
            fileArr->Add(digiFile);
            treeArr->Add(digiTree);
        } else {
            digiFile = 0;
            std::cout << "digi data file could not be opened!!" << std::endl;
        }
    }
    
    if (reconFile) {
        delete rec; 
        rec = 0;
        reconTree = 0;
        delete reconFile;
        reconFile = 0;
    }
    
    if (reconFileName != "") {
        reconFile = new TFile(reconFileName);
        if (reconFile->IsOpen() == kTRUE) {
            reconTree = (TTree*)gDirectory->Get("Recon");
            rec = 0;
            reconTree->SetBranchAddress("ReconEvent",&rec);
            fileArr->Add(reconFile);
            treeArr->Add(reconTree);
        } else {
            reconFile = 0;
            std::cout << "recon data file could not be opened!!" << std::endl;
        }
    }
    
    if (anaTupFile) {
        delete anaTup; 
        anaTupTree = 0;
        delete anaTupFile; 
        anaTupFile = 0;
    }
    
    if (anaTupFileName != "") {
        anaTupFile = new TFile(anaTupFileName);
        if (anaTupFile->IsOpen() == kTRUE) {
            anaTupTree = (TTree*)gDirectory->Get("1");
	    TObjArray* bList = anaTupTree->GetListOfBranches();
	    for (Long64_t ib=0; ib<bList->GetEntries(); ib++) {
	      const char* name = bList->At(ib)->GetName();
	      anaTupTree->SetBranchAddress(name,&anaTup[ib]);
	    }
            fileArr->Add(anaTupFile);
            treeArr->Add(anaTupTree);
        } else {
            anaTupFile = 0;
            std::cout << "anaTup data file could not be opened!!" << std::endl;
        }
    }

    m_StartEvent = 0;
    pruned = 0;
}


inline UInt_t RootTreeAnalysis::GetEvent(Long64_t ievt) {
    // Purpose and Method:  Get the event, ievt, for all trees
    //    We could be processing single files or chains, 
	//    This routine handles both casees.

    // if using regular trees - we check the array of open trees and
    // move the event pointer to the requested event
    UInt_t nb = 0;
    if (treeArr) {
        for (Long64_t i = 0; i < treeArr->GetEntries(); i++) {
            nb += ((TTree*)treeArr->At(i))->GetEvent(ievt);
        }
        return nb;
    }
    
    // if using chains, check the array of chains and move
    // the event pointer to the requested event
    if (chainArr) {
        for (Long64_t i = 0; i < chainArr->GetEntries(); i++) {
            nb += ((TChain*)chainArr->At(i))->GetEvent(ievt);
        }
        return nb;
    }

  return nb;
}



inline Long64_t RootTreeAnalysis::GetEntries() const {    
    // Purpose and Method:  Determine the number of events to iterate over
    //   checking to be sure that the requested number of events is less than
    //   the min number of events in all files

    Long64_t nentries = 0;
    if (treeArr) {
        nentries = ((TTree*)treeArr->At(0))->GetEntries();
        for (Long64_t i = 1; i < treeArr->GetEntries(); i++) {
            nentries = TMath::Min(nentries, ((TTree*)treeArr->At(i))->GetEntries());
        }
        return nentries;
    }
    
    if (chainArr) {
        nentries = ((TChain*)chainArr->At(0))->GetEntries();
        for (Long64_t i = 1; i < chainArr->GetEntries(); i++) {
            nentries = TMath::Min(nentries, ((TChain*)chainArr->At(i))->GetEntries());
        }
        return nentries;
    }
    
    return nentries;
}


inline void RootTreeAnalysis::MakeHistList() {
    // Purpose and Method:  Make a THashList of histograms
    //   This avoids the need to refresh the histogram pointers
    
    if (m_histList) delete m_histList;
    
    m_histList = new THashList(30, 5);
    
    TList* list = histFile->GetList();
    TIter iter(list);
    
    TObject* obj = 0;
    
    while (obj=iter.Next()) {
        m_histList->Add(obj);
    }
}

inline void RootTreeAnalysis::HistClear() {
    // Purpose and Method:  Clear histograms by iterating over the THashList
    
    if (!m_histList) return;
    
    TIter iter(m_histList);
    
    TObject* obj = 0;
    
    while ( obj=(TObject*)iter.Next() ) {
        ((TH1*)obj)->Reset();        
    }
}

inline void RootTreeAnalysis::Clear() {
    histFile = 0; 
    m_histList = 0;
    
    digiFile = 0; 
    reconFile = 0;
    mcFile = 0;
    anaTupFile = 0;
    
    digiTree = 0; 
    reconTree = 0;
    mcTree = 0;
    anaTupTree = 0;
    
    m_digiChain = 0;
    m_recChain = 0;
    m_mcChain = 0;
    m_anaTupChain = 0;
    
    evt = 0;
    rec = 0;
    mc = 0;
    for (int i=0; i<1000; i++) {
      anaTup[i]=0;
    }
    
    fileArr = 0;
    treeArr = 0;
    chainArr = 0;
}

#endif

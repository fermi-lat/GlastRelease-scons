/** @class TreeMaker
 * @brief This class reads a Digi, and a Recon Root file and merges their
 * contents into a single TTree.
 *
 * It is meant for single tower analysis.
 *
 */

#ifndef TreeMaker_h
#define TreeMaker_h 1

#if !defined(__CINT__)
//Need these includes if we wish to compile this code
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

class TreeMaker {
public :

    TFile       *TreeFile;
    /// Input digitization file
    TFile       *digiFile;   
    /// Input reconstruction file
    TFile       *reconFile;  
    /// Input monte carlo file
    TFile       *mcFile;     
    /// pointer to the digi tree
    TTree       *digiTree;    
    /// pointer to the reconstruction tree
    TTree       *reconTree;
    /// pointer to the monte carlo tree
    TTree       *mcTree;      
    /// Optional TChain input
    TChain      *m_digiChain, *m_recChain, *m_mcChain;
    /// pointer to a DigiEvent
    DigiEvent   *evt;
    /// pointer to a ReconEvent
    ReconEvent  *rec;
    /// Pointer to a McEvent
    McEvent     *mc;
    /// name of the output Tree ROOT file
    char        *m_TreeFileName; 
    /// Arrays that contain pointers to the TFile, TTree, and TChains
    TObjArray   *fileArr, *treeArr, *chainArr;

    TObjArray *TreeCollection;
    /// Default ctor, requires that that user calls TreeMaker::Init
    /// to setup access to specific ROOT files.
    TreeMaker(); 
        
	/// Standard ctor, where user provides the names of the input root files
	/// and optionally the name of the output ROOT Tree file
    TreeMaker( 
        const char* digiFileName, 
        const char* reconFileName="", 
        const char* mcFileName="", 
        const char* TreeFileName=""); 

	/// Special ctor which accepts TChains for input files
    TreeMaker( 
        TChain *digiChain, 
        TChain *recChain = 0, 
        TChain *mcChain = 0, 
        char* TreeFileName="MyRootFile.root");

    ~TreeMaker();  

    /// allow the user to specify their own file name for the output ROOT file
    void SetTreeFileName(const char *TreeFileName) { m_TreeFileName = const_cast<char*>(TreeFileName); }
    /// re-init with these input Root files
    void Init(  const char* digiFileName="", 
		const char* reconFileName="", 
		const char* mcFileName=""); 
    
    /// process events
    void CreateTree(Int_t numEvents=100000); 
    void CreateDigiTree(Int_t numEvents=100000); 
    
    /// returns number of events in all open files
    UInt_t GetEntries() const;
    /// retrieve a pointer to event
    UInt_t GetEvent(UInt_t, UInt_t, UInt_t);
    
private:
    /// reset all member variables
    void Clear(); 

    /// event processing for the monte carlo data
    void McData();

};


inline TreeMaker::TreeMaker() 
{
    Clear();
}

inline TreeMaker::TreeMaker(const char* digiFileName, 
			    const char* reconFileName, 
			    const char* mcFileName, 
			    const char* TreeFileName)
{
  // Purpose and Method:  Standard constructor where the user provides the 
  //  names of input ROOT files and optionally the name of the output ROOT
  // file
    printf(" opening files:\n\tdigi:\t%s\n\trecon:\t%s\n\tmc:\t%s\n",
	   digiFileName, reconFileName, mcFileName);
    
    Clear();
    
    SetTreeFileName(TreeFileName);
    
    Init(digiFileName, reconFileName, mcFileName);
    
}

inline TreeMaker::TreeMaker(TChain *digiChain, 
			    TChain *recChain, 
			    TChain *mcChain, 
			    char* TreeFileName)
{
    Clear();
    
    SetTreeFileName(TreeFileName);
    
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
    
}


inline TreeMaker::~TreeMaker() {
  //    TreeFile->Close();
    
    if (TreeFile) delete TreeFile;
    
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
    
    if (fileArr) delete fileArr;
    if (treeArr) delete treeArr;
    if (chainArr) delete chainArr;

	Clear();
}



inline void TreeMaker::Init(const char* digiFileName, const char* reconFileName, const char* mcFileName)
{
    // Purpose and Method:  Re-initialize file, tree, event pointers, using the 
    
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

    if ( strlen(mcFileName) ) {
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

    if ( strlen(digiFileName) ) {
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
    
    if ( strlen(reconFileName) ) {
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
    //    std::cout<<"---------** evt = "<<evt<<" ,rec = "<<rec<<" ,mc = "<<mc<<std::endl;
}


inline UInt_t TreeMaker::GetEvent(UInt_t mc, UInt_t digi, UInt_t recon) {
    // Purpose and Method:  Get the event, iEvt, for all trees
    //    We could be processing single files or chains, 
    //    This routine handles both casees.

    // if using regular trees - we check the array of open trees and
    // move the event pointer to the requested event
    UInt_t nb = 0;
    UInt_t iEvt = 0;
    if (treeArr) {
        for ( Int_t i=0; i<treeArr->GetEntries(); i++ ) {
            TTree* t = (TTree*)treeArr->At(i);
            // I am sure this could be done more elegant!
            TString name = t->GetName();
            if ( name == "Digi" )
                iEvt = digi;
            else if ( name == "Recon" )
                iEvt = recon;
            else if ( name == "Mc" )
                iEvt = mc;
            else {
                std::cerr << "TreeMaker::GetEvent: don't know what to do with tree " << name << std::endl;
                std::exit(42);
            }
            nb += ((TTree*)treeArr->At(i))->GetEvent(iEvt);
        }
        return nb;
    }

    // the code for chains has to get modified too!
    std::cerr << "TreeMaker::GetEvent: don't know what to do with chains!" << std::endl;
    std::exit(42);
    
    // if using chains, check the array of chains and move
    // the event pointer to the requested event
    if (chainArr) {
        for (Int_t i = 0; i < chainArr->GetEntries(); i++) {
            nb += ((TChain*)chainArr->At(i))->GetEvent(iEvt);
        }
        return nb;
    }

  return nb;
}


inline UInt_t TreeMaker::GetEntries() const {    
    // Purpose and Method:  Determine the number of events to iterate over
    //   checking to be sure that the requested number of events is less than
    //   the min number of events in all files

    UInt_t nentries = 0;
    if (treeArr) {
        nentries = (UInt_t)((TTree*)treeArr->At(0))->GetEntries();
        for (Int_t i = 1; i < treeArr->GetEntries(); i++) {
            nentries = TMath::Min(nentries, (UInt_t)((TTree*)treeArr->At(i))->GetEntries());
        }
        return nentries;
    }
    
    if (chainArr) {
        nentries = (UInt_t)((TChain*)chainArr->At(0))->GetEntries();
        for (Int_t i = 1; i < chainArr->GetEntries(); i++) {
            nentries = TMath::Min(nentries, (UInt_t)((TChain*)chainArr->At(i))->GetEntries());
        }
        return nentries;
    }
    
    return nentries;
}


inline void TreeMaker::Clear() {
    TreeFile = 0; 
    
    digiFile = 0; 
    reconFile = 0;
    mcFile = 0;
    
    digiTree = 0; 
    reconTree = 0;
    mcTree = 0;
    
    m_digiChain = 0;
    m_recChain = 0;
    m_mcChain = 0;
    
    evt = 0;
    rec = 0;
    mc = 0;
    
    fileArr = 0;
    treeArr = 0;
    chainArr = 0;
}

#endif

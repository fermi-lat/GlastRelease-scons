///////////////////////////////////////////////////////////
/*!
\class RootTreeAnalysis
\brief   This class can manipulate a digi and a recon Root file
at the same time.  This class is intended to provide useful manipulation
of the Root Event loop. Users put their analysis code
into the Go function (in RootTreeAnalysis.cxx). They should not
need to look at RootTreeAnalysis.h (except to see the interface)

allows for:
init and re-init use of a Root file
clear all histograms
'go n events' allowing to continue on in the file
or 'rewind'

  Example of use:

  gROOT->LoadMacro("RootTreeAnalysis.cxx");     // 'compile' class
  RootTreeAnalysis m("myDigiFile.root", "myReconFile.root"); // create RootTreeAnalysis object
  m.Go(500);      // loop over 500 events. Go contains your analysis code
  ... look at histograms ...
  m.Go()          // look at remainder of file
  ... look at histograms ...
  m.HistClear();      // clear histograms
  m.Init("AnotherRootFile.root");
  m.Go(50);
  ... and so on ...

After editing your Go function, you need to issue a gROOT->Reset() and
repeat the above sequence starting from the .L RootTreeAnalysis.cxx.

If you only have a digi or only a recon root file... setup RootTreeAnalysis like this:
if you only have a recon root file:
RootTreeAnalysis *m = new RootTreeAnalysis("", "myReconFile.root")
if you only have a digi root file:
RootTreeAnalysis *m = new RootTreeAnalysis("myDigiFile.root", "")

Version 0.1 17-Mar-1999 Richard Creation
Version 1.0 Spring, 2000 Revised for use with GLAST 1999 TestBeam
Version 1.5 25-Mar-2001  Revised for use with GLAST 2001 Balloon
Version 2.0 14-Aug-2001 Final version for GLAST 2001 Balloon
*/


#ifndef RootTreeAnalysis_h
#define RootTreeAnalysis_h

#include "iostream.h"

// Forward declarations
class DigiEvent;
class Recon;
class McEvent;

class RootTreeAnalysis {
public :
    TFile       *histFile;    // histogram file
    TFile       *digiFile;     // input digi file
    TFile       *reconFile;   // input recon file
    TFile       *mcFile;      // input mc file
    TTree       *digiTree;     //pointer to the digi TTree
    TTree       *reconTree;   //pointer to the recon TTree
    TTree       *mcTree;      
    TChain      *m_digiChain, *m_recChain, *m_mcChain;
    DigiEvent   *evt;
    ReconEvent       *rec;
    McEvent     *mc; // disable until MC Root classes are sorted out
    char        *m_histFileName; // actual name of the output ROOT file

                
    TObjArray   *fileArr, *treeArr, *chainArr;

    /// default constructor
    RootTreeAnalysis(); 
    /// ctr with root file name
    RootTreeAnalysis(   const char* digiFileName="", 
                        const char* reconFileName="", 
                        const char* mcFileName="", const char *histtFileName="Histograms.root"); 
    RootTreeAnalysis( TChain *digiChain = 0, TChain *recChain = 0, TChain *mcChain = 0, const char *histFileName="Histograms.root");
    ~RootTreeAnalysis();  
    /// start next Go with this event
    void StartWithEvent(Int_t event) { m_StartEvent = event; };  
    /// reset for next Go to start at beginning of file 
    void Rewind() { m_StartEvent = 0; }; 
    /// allow the user to specify their own file name for the output ROOT file
    void SetHistFileName(const char *histFileName) { m_histFileName = histFileName; }
    /// re-init with these root files
    void Init(  const char* digiFileName="", 
                const char* reconFileName="", 
                const char* mcFileName=""); 
    /// Initialize the pedestal variables to zero
    void ZeroPeds(); 
    /// read the pedestal files, and fill the local pedestal variables
    void ReadPedFiles();
    /// prints out the pedestals for CAL, ACD, and XGTs
    void PrintPeds();
    /// initialize the threshold array to zero
    void ZeroVetoThresh();
    /// read in the ACD veto thresholds for the tiles
    void ReadVetoThresh();
    /// prints out the current set of ACD veto thresholds
    void PrintVetoThresh();
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
    /// loop over events
    void Go(Int_t numEvents=100000); 
    /// returns min num of events in all open files
    UInt_t GetEntries();
    /// retrieve a pointer to event number.
    UInt_t GetEvent(UInt_t ievt);


private:
    /// starting event
    Int_t m_StartEvent;
    /// list of user histograms
    THashList *m_histList;

    // storage for pedestals
    Float_t m_calPed[4][8][10][2];
    Float_t m_acdPed[1001]; 
    Float_t m_xgtPed[12]; 
    // stores the ACD veto threshold
    Float_t m_acdVeto[1001];

    /// reset all member variables
    void Clear();  
    void DigiHistDefine();
    void ReconHistDefine();

    void McHistDefine();

    void DigiCalHistDefine();
    void DigiTkrHistDefine();
    void DigiAcdHistDefine();

    void ReconCalHistDefine();
    void ReconTkrHistDefine();
    void ReconAcdHistDefine();

    /// event processing for the digi data
    void DigiTkr();
    void DigiCal();
    void DigiAcd();

    /// event processing for the reconstruction data
    void ReconTkr();
    void ReconCal();
    void ReconAcd();
        
    /// event processing for the monte carlo data
    void McData();
};

#endif

RootTreeAnalysis::~RootTreeAnalysis() {
    // destructor: cleanup

    histFile->Close();

    //if (m_histList) delete m_histList;

    if (histFile) delete histFile;

    if (digiFile) delete digiFile;
    if (reconFile) delete reconFile;
    if (mcFile) delete mcFile;

    if (evt) delete evt;
    if (rec) delete rec;
    if (mc) delete mc;

    digiTree = 0;
    reconTree = 0;
    mcTree = 0;

    if (fileArr) delete fileArr;
    if (treeArr) delete treeArr;
    if (chainArr) delete chainArr;

}

void RootTreeAnalysis::Clear() {
    histFile=0; 
    m_histList = 0;

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
    ZeroPeds();
    ZeroVetoThresh();
}

RootTreeAnalysis::RootTreeAnalysis() 
{
    Clear();
    ReadPedFiles();
    ReadVetoThresh();
}

RootTreeAnalysis::RootTreeAnalysis(const char* digiFileName, 
                                   const char* reconFileName, 
                                   const char* mcFileName, 
                                   const char* histFileName)
{
    // constructor: set up root file and init things

    printf(" opening files:\n\tdigi:\t%s\n\trecon:\t%s\n\tmc:\t%s\n",digiFileName, reconFileName, mcFileName);

    Clear();

    SetHistFileName(histFileName);
    HistDefine();
    MakeHistList();

    Init(digiFileName, reconFileName, mcFileName);

//    ReadPedFiles();
//   ReadVetoThresh();
}

RootTreeAnalysis::RootTreeAnalysis(TChain *digiChain, 
                                   TChain *recChain, 
                                   TChain *mcChain, 
                                   const char *histFileName) {
    Clear();

    SetHistFileName(histFileName);
    HistDefine();
    MakeHistList();

    if (chainArr) delete chainArr;
    chainArr = new TObjArray();

    if (digiChain != 0) {
        evt = new DigiEvent();
        m_digiChain = digiChain;
        m_digiChain->SetBranchAddress("DigiEvent",&evt);
        chainArr->Add(m_digiChain);
    }

    if (recChain != 0) {
        m_recChain = recChain;
        rec = new Recon();
        m_recChain->SetBranchAddress("Recon",&rec);
        chainArr->Add(m_recChain);
    }

    if (mcChain != 0) {
        m_mcChain = mcChain;
        mc = new McEvent();
        m_mcChain->SetBranchAddres("mc",&mc);
        chainArr->Add(m_mcChain);
    }

    ReadPedFiles();
    ReadVetoThresh();
}



void RootTreeAnalysis::Init(const char* digiFileName, const char* reconFileName, const char* mcFileName)
{
// re-initialize file, tree, event. Histograms are *not* cleared.

    if (fileArr) delete fileArr;
    fileArr = new TObjArray();

    if (treeArr) delete treeArr;
    treeArr = new TObjArray();

//   Set branch addresses

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
//            mc = new McEvent();
            mc = 0;
            mcTree->SetBranchAddress("McEvent",&mc);
            fileArr->Add(mcFile);
            treeArr->Add(mcTree);
        } else {
            mcFile = 0;
            std::cout << "mc data file could not be opened!!" << std::endl;
        }
    }

    m_StartEvent = 0;

}


UInt_t RootTreeAnalysis::GetEvent(UInt_t ievt) {

    // Either using a single file and tree or a chain of files

    // if using regular trees - we check the array of open trees and
    // move the event pointer to the requested event
    UInt_t nb = 0;
    if (treeArr) {
        for (Int_t i = 0; i < treeArr->GetEntries(); i++) {
            nb += ((TTree*)treeArr->At(i))->GetEvent(ievt));
        }
        return nb;
    }

    // if using chains, check the array of chains and move
    // the event pointer to the requested event
    if (chainArr) {
        for (Int_t i = 0; i < chainArr->GetEntries(); i++) {
            nb += ((TChain*)chainArr->At(i))->GetEvent(ievt));
        }
        return nb;
    }

}


UInt_t RootTreeAnalysis::GetEntries() {

    UInt_t nentries = 0;
    // Determine the number of events to iterate over - checking to be sure that
    // the requested number of events is less than the min number of events in all files
    if (treeArr) {
        nentries = ((TTree*)treeArr->At(0))->GetEntries();
        for (Int_t i = 1; i < treeArr->GetEntries(); i++) {
            nentries = TMath::Min(nentries, ((TTree*)treeArr->At(i))->GetEntries());
        }
        return nentries;
    }

    if (chainArr) {
        nentries = ((TChain*)chainArr->At(0))->GetEntries();
        for (Int_t i = 1; i < chainArr->GetEntries(); i++) {
            nentries = TMath::Min(nentries, ((TChain*)chainArr->At(i))->GetEntries());
        }
        return nentries;
    }

}


void RootTreeAnalysis::MakeHistList() {
// make a THashList of histograms

     if (m_histList) delete m_histList;

     TFile *histFile = (TFile*)gROOT->GetFile(m_histFileName);
     m_histList = new THashList(30, 5);

     TList* list = histFile->GetList();
     TIter*  iter = new TIter(list);

     TObject* obj = 0;

     while (obj=iter->Next()) {
       
       //  if (obj->InheritsFrom("TH1")) {
             m_histList->Add(obj);
       //  }
     }
     delete iter;
}

void RootTreeAnalysis::HistClear() {
// clear histograms

    if (!m_histList) return;

    TIter iter(m_histList);
    
    TObject* obj = 0;
    
    while ( obj=(TObject*)iter.Next() ) {
        ((TH1*)obj)->Reset();        
    }


}


void RootTreeAnalysis::ZeroPeds() {

    for (int i = 0; i < 4; i++) {
        for (int lyr = 0; lyr < 8; lyr++){
            for (int col = 0; col < 10; col++) {
                for (int side = 0; side < 2; side++) {
                    m_calPed[i][lyr][col][side]=0.0;
                }
            }
        }
    }

    for (int i=0;i<1001;i++) m_acdPed[i]=0.0; 
    for (int i=0;i<12;i++) m_xgtPed[i]=0.0; 

}

void RootTreeAnalysis::ZeroVetoThresh() {
    for (int i=0; i<1001;i++) m_acdVeto[i] = 0.0;
}

void RootTreeAnalysis::ReadPedFiles() {
    // Read in the pedestal files and store the data in arrays
    
    //start reading CAL pedestal file
    int ncols,nlines=0;
    int facd,flyr,fcol,fside;
    float fped,fsig;
    cout<<"Read in CAL pedestal file"<<endl;        
    char *RootAnalysisPath = ::getenv("ROOTANALYSIS");
    char fileName[800];
    sprintf(fileName, "%s/calibration/bfem_peds_cal.dat", RootAnalysisPath);
    FILE *fCal = fopen(fileName,"r");
    
    while (1)
    {
        ncols = fscanf(fCal,"%d %d %d %d %f %f",&facd,&flyr,&fcol,&fside,&fped,&fsig);
        if (ncols<0) break;
        m_calPed[facd][flyr][fcol][fside]=fped;
        ++nlines;
    }
    
    fclose(fCal);
    cout<<"Done"<<endl;        
    
    //start reading ACD ped file
    cout<<"Read in ACD pedestal file"<<endl; 
    
    ncols=0; 
    nlines=0; 
    sprintf(fileName, "%s/calibration/bfem_peds_acd.dat", RootAnalysisPath);
    FILE *fAcd = fopen(fileName,"r");  
    while (1)  
    {  
        ncols = fscanf(fAcd,"%d %f",&facd,&fped);  
        if (ncols<0) break;  
        m_acdPed[facd]=fped;  
        ++nlines;  
    }  
    
    fclose(fAcd);  
    cout<<"Done"<<endl; 
    
    //start reading XGT ped file
    cout<<"Read in XGT pedestal file"<<endl; 
    
    ncols=0; 
    nlines=0; 
    sprintf(fileName, "%s/calibration/bfem_peds_xgt.dat", RootAnalysisPath);
    FILE *fXgt = fopen(fileName,"r");  
    while (1)  
    {  
        ncols = fscanf(fXgt,"%d %f",&facd,&fped);  
        if (ncols<0) break;  
        m_xgtPed[facd]=fped;  
        ++nlines;  
    }  
    
    fclose(fXgt);  
    cout<<"Done"<<endl; 
}

void RootTreeAnalysis::ReadVetoThresh() {
    char *RootAnalysisPath = ::getenv("ROOTANALYSIS");
    char fileName[800];
    sprintf(fileName, "%s/calibration/bfem_vetothresh_acd.dat", RootAnalysisPath);
    FILE *fAcd = fopen(fileName,"r");

    int counter;
    Int_t index;
    Float_t thresh;
    for (counter = 0; counter < 13; counter++) {
        fscanf(fAcd, "%d %f", &index, &thresh);
        m_acdVeto[index] = thresh;
    }

    fclose(fAcd);
}

void RootTreeAnalysis::PrintPeds() {
    cout << "The Pedestal values are: " << endl;
    cout << "CAL" << endl;

    cout << "ACD Id\tPedestal" << endl;
    cout << "0\t" << m_acdPed[0] << endl;
    cout << "1\t" << m_acdPed[1] << endl;
    cout << "10\t" << m_acdPed[10] << endl;
    cout << "11\t" << m_acdPed[11] << endl;
    cout << "100\t" << m_acdPed[100] << endl;
    cout << "110\t" << m_acdPed[110] << endl;
    cout << "200\t" << m_acdPed[200] << endl;
    cout << "210\t" << m_acdPed[210] << endl;
    cout << "300\t" << m_acdPed[300] << endl;
    cout << "310\t" << m_acdPed[310] << endl;
    cout << "400\t" << m_acdPed[400] << endl;
    cout << "410\t" << m_acdPed[410] << endl;
    cout << "1000\t" << m_acdPed[1000] << endl;

    cout << "XGT" << endl;
    cout << "2000\t" << m_xgtPed[0] << endl;
    cout << "2001\t" << m_xgtPed[1] << endl;
    cout << "2010\t" << m_xgtPed[10] << endl;
    cout << "2011\t" << m_xgtPed[11] << endl;

}

void RootTreeAnalysis::PrintVetoThresh() {

    cout << endl << "The ACD Veto Thresholds are: " << endl;
    cout << "ACD Id\tVeto Threshold" << endl;
    cout << "0\t" << m_acdVeto[0] << endl;
    cout << "1\t" << m_acdVeto[1] << endl;
    cout << "10\t" << m_acdVeto[10] << endl;
    cout << "11\t" << m_acdVeto[11] << endl;
    cout << "100\t" << m_acdVeto[100] << endl;
    cout << "110\t" << m_acdVeto[110] << endl;
    cout << "200\t" << m_acdVeto[200] << endl;
    cout << "210\t" << m_acdVeto[210] << endl;
    cout << "300\t" << m_acdVeto[300] << endl;
    cout << "310\t" << m_acdVeto[310] << endl;
    cout << "400\t" << m_acdVeto[400] << endl;
    cout << "410\t" << m_acdVeto[410] << endl;
    cout << "1000\t" << m_acdVeto[1000] << endl;
}

#ifndef __LEANINGTOWER_EVENT__
#define __LEANINGTOWER_EVENT__

#include "Recon.h"
#include "Layer.h"

#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TGraph.h"
#include "TString.h"
#include "TObjString.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class Event {
 public:
    Event(TString filename, TList* geometry);
    virtual ~Event();

    void Go(int event) {
        SelectedEvent = event;
        myTree->GetEntry(SelectedEvent);
    }

    int GetEntries()           const { return (int)myTree->GetEntries(); }
    Int_t GetEventId()         const { return EventId; }
    Int_t GetRunId()           const { return RunId; }
    Int_t GetTkrTotalNumHits() const { return TkrTotalNumHits; }
    Double_t GetEbfTime()      const { return EbfTime; }

    int GetPlaneNumHits(TString PlaneName);
    int *GetPlaneHits(TString PlaneName);
    Bool_t GetTriggerReq(TString PlaneName, Bool_t side);

    // hits and clusters
    TGraph GetTGraphHits(TString view);
    TGraph GetTGraphHits(int view);
    std::vector<double> GetClusters(TString PlaneName);
    TGraph GetTGraphClusters(TString view);
    TGraph GetTGraphClusters(int view);

    UInt_t GetTkr3RowBits() { return Tkr3RowBits;}

    // Recon

    Recon* GetRecon() {
        myRecon->GetEvent(SelectedEvent);
        return myRecon;
    }

 private:
    int NumberOfEvents;
    static const int GAP = 10;
    int SelectedEvent;
    Int_t TkrTotalNumHits, RunId, EventId;
    Double_t EbfTime;
    int* TkrHits;

    TTree* myTree;
    TFile* myFile;
    TList*  myGeometry;
    Recon* myRecon;

    UInt_t Tkr3RowBits;
  
    // last line of class def
    ClassDef(Event, 1)
};

#endif

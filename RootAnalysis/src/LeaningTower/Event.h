#include "TFile.h"
#include "TTree.h"
#include "TMap.h"
#include "TGraph.h"
#include "TString.h"
#include "TObjString.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>


class Event
{
 public:
  Event(TString filename, TMap *geometry);
  inline void Go(int event)
    {
      SelectedEvent=event;
    }
  
  inline  int GetTkrTotalNumHits() 
    {
      myTree->GetEntry(SelectedEvent);
      return TkrTotalNumHits;
    }

  inline  int GetRunId() 
    {
      myTree->GetEntry(SelectedEvent);
      return RunId;
    }
  
  inline  int GetEventId()
    {
      myTree->GetEntry(SelectedEvent);
      return EventId;
    }
  Double_t GetEbfTime() { return EbfTime; }
  
  int GetLayerNumHits(TString LayerName);
  int *GetLayerHits(TString LayerName);
  Bool_t Event::GetTriggerReq(TString LayerName, Bool_t side);
  std::vector<double> GetClusters(TString LayerName);
  
 private:
  int NumberOfEvents;
  static const int GAP=10;
  int SelectedEvent;
  int TkrTotalNumHits, RunId, EventId;
  Double_t EbfTime;
  //  int TkrNumHits;
  int *TkrHits;
  

  TTree *myTree;
  TFile *myFile;
  TMap * myGeometry;
  
  
};

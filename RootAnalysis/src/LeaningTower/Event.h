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
  
  int GetLayerNumHits(TString LayerName);
  int *GetLayerHits(TString LayerName);
  std::vector<double> GetClusters(TString LayerName);
  
 private:
  int NumberOfEvents;

  int SelectedEvent;
  int TkrTotalNumHits, RunId, EventId;
  //  int TkrNumHits;
  int *TkrHits;

  TTree *myTree;
  TFile *myFile;
  TMap * myGeometry;
  
  
};

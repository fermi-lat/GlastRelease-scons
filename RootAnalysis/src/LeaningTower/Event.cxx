#include "Event.h"

#include "Layer.h"
#include "Recon.h"

#include <vector>

Event::Event(TString filename, TMap* geometry) {
    myFile = new TFile(filename);
    myGeometry = geometry;  

    myTree = (TTree*)myFile->Get("Header");
    myTree->SetBranchAddress("EventId",&EventId);
    myTree->SetBranchAddress("RunId",&RunId);
    myTree->SetBranchAddress("TkrTotalNumHits",&TkrTotalNumHits);
    myTree->SetBranchAddress("EbfTime",&EbfTime);

    NumberOfEvents = (int)myTree->GetEntries();
    std::cout << "Number of Events: " << NumberOfEvents << std::endl;
  
    TMapIter ti(myGeometry);
    TObjString* key;
    while ( ( key = (TObjString*)ti.Next() ) ) {
      Layer *aLayer = (Layer*)myGeometry->GetValue(key);
      aLayer->SetTree(myFile);
    }

    myRecon = new Recon(myFile);
}

Event::~Event() {
    delete myRecon;
    myFile->Close();
    delete myFile;
}

int Event::GetLayerNumHits(TString LayerName)
{
  Layer *aLayer = ((Layer*) myGeometry->GetValue(LayerName));
  aLayer->GetEvent(SelectedEvent);
  return aLayer->TkrNumHits;
}

int *Event::GetLayerHits(TString LayerName)
{
  Layer *aLayer = ((Layer*) myGeometry->GetValue(LayerName));
  aLayer->GetEvent(SelectedEvent);
  return aLayer->TkrHits;
}

Bool_t Event::GetTriggerReq(TString LayerName, Bool_t side) {
  Layer *aLayer = (Layer*)myGeometry->GetValue(LayerName);
  aLayer->GetEvent(SelectedEvent);
  return aLayer->GetTriggerReq(side);
}    

TGraph Event::GetTGraphHits(TString view) {
    if ( view == "X" )
        return GetTGraphHits(0);
    else if ( view == "Y" )
        return GetTGraphHits(1);
    else
        std::cerr << "Event::GetTGraphHits: view = " << view << std::endl;

    std::exit(42);
    return TGraph();
}

TGraph Event::GetTGraphHits(int view) {
    TGraph tg;
    TMapIter ti(myGeometry);
    TObjString* key;

    while ( ( key = (TObjString*)ti.Next() ) ) {
        Layer* l = (Layer*)myGeometry->GetValue(key);
        if ( view == 0 && !l->IsX() )
            continue;
        if ( view == 1 && !l->IsY() )
            continue;

        TString layerName = key->String();
        int layerNumHits = GetLayerNumHits(layerName);
        if ( layerNumHits <= 0 )
            continue;

        double z = l->GetHeight();
        int* layerHits = GetLayerHits(layerName);
      
        for ( int i=0; i<layerNumHits; ++i ) {
            double pos = l->GetCoordinate(layerHits[i]);
            //            std::cout << layerName << ' ' << view << ' ' << pos << ' ' << z << std::endl;
            tg.SetPoint(tg.GetN(), pos, z);
        }
    }
    return tg;
}

std::vector<double> Event::GetClusters(TString LayerName)
{
  Layer *aLayer = ((Layer*) myGeometry->GetValue(LayerName));
  aLayer->GetEvent(SelectedEvent);
  
  int NumHits = aLayer->TkrNumHits;
  int *Hits   = aLayer->TkrHits;
  
  std::vector<double> Cluster;
  if(NumHits==0) return Cluster;
  //  Cluster->Set(NumHits);

  int LastClusterStrip = Hits[0];
  
  Cluster.push_back(aLayer->GetCoordinate(LastClusterStrip));
  
  int ClusterSize=1;
  int Clusters=0;
  
  for(int i=1; i<NumHits; i++)
    {
      int nextHit =  Hits[i];
      if(nextHit < LastClusterStrip+GAP) 
	{
	  Cluster[Clusters] = (ClusterSize * Cluster[Clusters] + 
			       aLayer->GetCoordinate(nextHit))/
	    ((double)(ClusterSize+1));
	  ClusterSize++;
	}
      else
	{
	  Clusters++;
	  LastClusterStrip  = nextHit;
	  Cluster.push_back(aLayer->GetCoordinate(LastClusterStrip));
	}
    }
  //  Cluster.Set(Clusters);
  // This is to have no more than one cluster per layer:

  if(Cluster.size()>1) 
    Cluster.clear();
  return Cluster;
}

TGraph Event::GetTGraphClusters(TString view) {
    if ( view == "X" )
        return GetTGraphClusters(0);
    else if ( view == "Y" )
        return GetTGraphClusters(1);
    else
        std::cerr << "Event::GetTGraphClusters: view = " << view << std::endl;

    std::exit(42);
    return TGraph();
}

TGraph Event::GetTGraphClusters(int view) {
    TGraph tg;
    TMapIter ti(myGeometry);
    TObjString* key;

    while ( ( key = (TObjString*)ti.Next() ) ) {
        Layer* l = ((Layer*)myGeometry->GetValue(key));
        if ( view == 0 && !l->IsX() )
            continue;
        if ( view == 1 && !l->IsY() )
            continue;

        TString layerName = key->String();
        int layerNumHits = GetLayerNumHits(layerName);
        if ( layerNumHits <= 0 )
            continue;

        double z = l->GetHeight();
        // getting cluster(s) from each layer
        std::vector<double> clusterPos = GetClusters(layerName);
        int clusterNum = clusterPos.size();
        for ( int i=0; i<clusterNum; ++i )
            tg.SetPoint(tg.GetN(), clusterPos[i], z);
    }
    return tg;
}

ClassImp(Event)

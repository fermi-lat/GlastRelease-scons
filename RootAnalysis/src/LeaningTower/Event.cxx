#include <vector>


#include "Event.h"
#include "Layer.h"


Event::Event(TString filename, TMap *geometry)
{
  myFile = new TFile(filename);
  myTree = (TTree*) myFile->Get("Header");
  
  //////////////////////////////////////////////////
  myGeometry = geometry;  
  //////////////////////////////////////////////////

  myTree->SetBranchAddress("EventId",&EventId);
  myTree->SetBranchAddress("RunId",&RunId);
  myTree->SetBranchAddress("TkrTotalNumHits",&TkrTotalNumHits);
  
  NumberOfEvents = myTree->GetEntries();
  std::cout << "Number of Events: " << NumberOfEvents << std::endl;
  //////////////////////////////////////////////////
  
  TMapIter ti(myGeometry);
  TObjString *key;
  while (key = (TObjString*)ti.Next())
    {
      Layer *aLayer = ((Layer*) myGeometry->GetValue(key));
      aLayer->SetTree(myFile);
    }
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
      if(nextHit < LastClusterStrip+5) 
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







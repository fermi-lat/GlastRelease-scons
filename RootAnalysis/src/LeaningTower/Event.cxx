#include "Event.h"

#include "Layer.h"
#include "Recon.h"

#include <vector>

Event::Event(TString filename, TList* geometry) {
    myFile = new TFile(filename);
    myGeometry = geometry;  

    myTree = (TTree*)myFile->Get("Header");
    myTree->SetBranchAddress("EventId",&EventId);
    myTree->SetBranchAddress("RunId",&RunId);
    myTree->SetBranchAddress("TemId",&TemId);
    myTree->SetBranchAddress("TkrTotalNumHits",&TkrTotalNumHits);
    myTree->SetBranchAddress("EbfTime",&EbfTime);
    myTree->SetBranchAddress("TkrDigi3RowBits",&TkrDigi3RowBits);
    myTree->SetBranchAddress("TkrTrgReq3RowBits",&TkrTrgReq3RowBits);

    NumberOfEvents = (int)myTree->GetEntries();
    std::cout << "Number of events: " << NumberOfEvents << std::endl;
  
    TIter next(myGeometry);
    while ( Layer* aPlane = (Layer*)next() )
        aPlane->SetTree(myFile);
    myRecon = new Recon(myFile);
}

Event::~Event() {
    delete myRecon;
    myFile->Close();
    delete myFile;
}

int Event::GetPlaneNumHits(TString PlaneName) {
    Layer* aPlane = (Layer*)myGeometry->FindObject(PlaneName);
    aPlane->GetEvent(SelectedEvent);
    return aPlane->TkrNumHits;
}

int* Event::GetPlaneHits(TString PlaneName) {
    Layer* aPlane = (Layer*)myGeometry->FindObject(PlaneName);
    aPlane->GetEvent(SelectedEvent);
    return aPlane->TkrHits;
}

Bool_t Event::GetTriggerReq(TString PlaneName, Bool_t side) {
    Layer* aPlane = (Layer*)myGeometry->FindObject(PlaneName);
    aPlane->GetEvent(SelectedEvent);
    return aPlane->GetTriggerReq(side);
}

int Event::GetToT(TString PlaneName, Bool_t side) {
    Layer* aPlane = (Layer*)myGeometry->FindObject(PlaneName);
    aPlane->GetEvent(SelectedEvent);
    return !side ? aPlane->ToT0 : aPlane->ToT1;
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
    TIter next(myGeometry);
    while ( Layer* aPlane = (Layer*)next() ) {
        if ( view == 0 && !aPlane->IsX() )
            continue;
        if ( view == 1 && !aPlane->IsY() )
            continue;

        const TString planeName = aPlane->GetName();
        const int planeNumHits = GetPlaneNumHits(planeName);
        if ( planeNumHits <= 0 )
            continue;

        int* planeHits = GetPlaneHits(planeName);
        for ( int i=0; i<planeNumHits; ++i )
            tg.SetPoint(tg.GetN(), aPlane->GetCoordinate(planeHits[i]),
                        aPlane->GetHeight());
    }
    return tg;
}

std::vector<float> Event::GetClusters(TString PlaneName) {
    Layer* aPlane = (Layer*)myGeometry->FindObject(PlaneName);
    aPlane->GetEvent(SelectedEvent);
  
    int NumHits = aPlane->TkrNumHits;
    int *Hits   = aPlane->TkrHits;
  
    std::vector<float> Cluster;
    if(NumHits==0) return Cluster;
    //  Cluster->Set(NumHits);

    int LastClusterStrip = Hits[0];
  
    Cluster.push_back(aPlane->GetCoordinate(LastClusterStrip));
  
    int ClusterSize=1;
    int Clusters=0;
  
    for ( int i=1; i<NumHits; i++ ) {
        int nextHit =  Hits[i];
        if ( nextHit < LastClusterStrip+GAP ) {
            Cluster[Clusters] = (ClusterSize * Cluster[Clusters] + 
                                 aPlane->GetCoordinate(nextHit))/
                (float)(ClusterSize+1);
            ClusterSize++;
	}
        else {
            Clusters++;
            LastClusterStrip  = nextHit;
            Cluster.push_back(aPlane->GetCoordinate(LastClusterStrip));
	}
    }
    //  Cluster.Set(Clusters);
    // This is to have no more than one cluster per plane:

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
    TIter next(myGeometry);
    while ( Layer* aPlane = (Layer*)next() ) {
        if ( view == 0 && !aPlane->IsX() )
            continue;
        if ( view == 1 && !aPlane->IsY() )
            continue;

        const TString planeName = aPlane->GetName();
        const int planeNumHits = GetPlaneNumHits(planeName);
        if ( planeNumHits <= 0 )
            continue;

        // getting cluster(s) from each plane
        std::vector<float> clusterPos = GetClusters(planeName);
        const int clusterNum = clusterPos.size();
        for ( int i=0; i<clusterNum; ++i )
            tg.SetPoint(tg.GetN(), clusterPos[i], aPlane->GetZ());
    }
    return tg;
}

ClassImp(Event)

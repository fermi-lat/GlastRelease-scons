#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TMap.h"
#include "TObjectTable.h"
#include "TObjString.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TText.h"
#include "TTree.h"

#include "Tracker.h"
#include "Layer.h"
#include "Event.h"
#include "Recon.h"

#include <algorithm>

#define TOWER 1

class EventDisplay {
private:
    Event*   myEvent;
    Tracker* myTracker;
    TCanvas* myEventDisplay;

    TText LabelNumHits[36];
    TText TriggerReqText[36][2];

    TGraph *TGraphXhits;
    TGraph *TGraphYhits;

    TGraph* TGraphXclusters;
    TGraph* TGraphYclusters;
    TGraph* TGraphXclusters1;
    TGraph* TGraphYclusters1;
    TGraph* TGraphXclusters2;
    TGraph* TGraphYclusters2;
    TLine anXtrack;
    TLine anYtrack;

public:
    EventDisplay(TString="MyRootFile.root");
    void Go(int=-1);
};


EventDisplay::EventDisplay(TString filename) {
    gStyle->SetCanvasColor(10);
    myTracker = new Tracker;
    if ( TOWER ) {  
        myTracker->loadGeometry(gSystem->ExpandPathName("$ROOTANALYSISROOT/src/LeaningTower/geometry/Tower0Geometry.txt"));
        myTracker->loadFitting(gSystem->ExpandPathName("$ROOTANALYSISROOT/src/LeaningTower/geometry/FittingPlanes.txt"));
        myTracker->IsTower(TOWER);
    }
    else
        myTracker->loadGeometry(gSystem->ExpandPathName("$ROOTANALYSISROOT/src/LeaningTower/geometry/stack2geometry.txt"));
  
    myEvent = new Event(filename, (TMap*)myTracker->GetGeometry());
  
    myEventDisplay = new TCanvas("myEventDisplay", "myEventDisplay", 1100, 800);  
    myTracker->Display(myEventDisplay);
  
    for ( int i=0; i<36; ++i ) {
        TriggerReqText[i][0].SetTextAlign(32);
        TriggerReqText[i][1].SetTextAlign(02);
        LabelNumHits[i].SetTextAlign(02);
        TriggerReqText[i][0].SetTextSize(0.02);
        TriggerReqText[i][1].SetTextSize(0.02);
        LabelNumHits[i].SetTextSize(0.02);
    }

    TGraphXhits = new TGraph;
    TGraphYhits = new TGraph;

    TGraphXclusters = new TGraph;
    TGraphYclusters = new TGraph;
    TGraphXclusters1 = new TGraph;
    TGraphYclusters1 = new TGraph;
    TGraphXclusters2 = new TGraph;
    TGraphYclusters2 = new TGraph;
}

void EventDisplay::Go(int numEvent) {
    // actually, numEvent is the number of the record in the tree, not the event id
    static int entry = -1;
    if ( numEvent < 0 ) // display the next (special case, the first)
        ++entry;
    else
        entry = numEvent;
    int lastEntry = myEvent->GetEntries() - 1;
    if ( entry > lastEntry ) {
        entry = lastEntry;
        std::cout << "displaying the last record of the tree" << std::endl;
    }

    myEvent->Go(entry);

    int TkrTotalNumHits = myEvent->GetTkrTotalNumHits();

    std::cout << "EventId         = " << myEvent->GetEventId() << std::endl;
    std::cout << "RunId           = " << myEvent->GetRunId() << std::endl;
    std::cout << "EbfTime         = " << myEvent->GetEbfTime() << std::endl;
    std::cout << "TkrTotalNumHits = " << TkrTotalNumHits << std::endl;

    Recon* recon = myEvent->GetRecon();
    int TkrNumClus = recon->GetTkrNumClus();
    int TkrNumTracks = recon->GetTkrNumTracks();
    int TkrTrk1NumClus = recon->GetTkrTrk1NumClus();
    std::cout << "TkrNumClus      = " << TkrNumClus << std::endl;
    std::cout << "TkrNumTracks    = " << TkrNumTracks << std::endl;
    std::cout << "TkrTrk1NumClus  = " << TkrTrk1NumClus << std::endl;
    if ( TkrTrk1NumClus > TkrNumClus || TkrNumClus > TkrTotalNumHits )
        std::cout << "??? TkrTrk1NumClus > TkrNumClus > TkrTotalNumHits ???" << std::endl;

    TMap *myGeometry = myTracker->GetGeometry();
  
    TMapIter ti(myGeometry);
    TObjString* key;

    int i=0;
    while ( ( key = (TObjString*)ti.Next() ) ) {
        Layer *aLayer = ((Layer*) myGeometry->GetValue(key));

        TString LayerName=key->String();

        double height = aLayer->GetHeight();
      
        int LayerNumHits = myEvent->GetLayerNumHits(LayerName);
        //      int *LayerHits   = myEvent->GetLayerHits(LayerName);
      
        //        TString title = (TString("(")+=LayerNumHits) + ")";

        if(aLayer->IsX()) 
            myEventDisplay->cd(1);      
        else
            myEventDisplay->cd(2);

        TriggerReqText[i][0].SetText( -5, height, aLayer->GetTriggerReq(false) ? "x" : ".");
        TriggerReqText[i][1].SetText(365, height, aLayer->GetTriggerReq(true) ? "x" : ".");
        LabelNumHits[i].SetText(410, height, (TString("(")+=LayerNumHits) + ")");
        TriggerReqText[i][0].Draw();
        TriggerReqText[i][1].Draw();
        LabelNumHits[i].Draw();
        i++;
    }

    delete TGraphXhits;
    TGraphXhits = new TGraph(myEvent->GetTGraphHits("X"));
    delete TGraphYhits;
    TGraphYhits = new TGraph(myEvent->GetTGraphHits("Y"));

    // TGraph of clusters of Nicola's clustering
    delete TGraphXclusters;
    TGraphXclusters = new TGraph(myEvent->GetTGraphClusters("X"));
    delete TGraphYclusters;
    TGraphYclusters = new TGraph(myEvent->GetTGraphClusters("Y"));

    // TGraph of all clusters found by recon
    delete TGraphXclusters1;
    TGraphXclusters1 = new TGraph(recon->GetAllClustersGraph(0));
    delete TGraphYclusters1;
    TGraphYclusters1 = new TGraph(recon->GetAllClustersGraph(1));

    // TGraph of all clusters on track1 found by recon
    delete TGraphXclusters2;
    TGraphXclusters2 = new TGraph(recon->GetTrk1ClustersGraph("X"));
    delete TGraphYclusters2;
    TGraphYclusters2 = new TGraph(recon->GetTrk1ClustersGraph("Y"));

    // Drawing
  
    TGraphXhits->SetMarkerStyle(7);
    TGraphYhits->SetMarkerStyle(7);
    TGraphXclusters->SetMarkerStyle(4); 
    TGraphYclusters->SetMarkerStyle(4);
    TGraphXclusters->SetMarkerColor(3); 
    TGraphYclusters->SetMarkerColor(3);
    TGraphXclusters1->SetMarkerStyle(5); 
    TGraphYclusters1->SetMarkerStyle(5);
    TGraphXclusters1->SetMarkerColor(3); 
    TGraphYclusters1->SetMarkerColor(3);
    TGraphXclusters2->SetMarkerStyle(4); 
    TGraphYclusters2->SetMarkerStyle(4);
    TGraphXclusters2->SetMarkerColor(3);
    TGraphYclusters2->SetMarkerColor(3);
  
    // display

    myEventDisplay->cd(1);
    TGraphXhits->Draw("P");
    /*
      if ( TGraphXclusters->GetN() > 0 )
      TGraphXclusters->Draw("P");
      anXtrack = Reconstruct(TGraphXclusters);
    */
    if ( TGraphXclusters1->GetN() > 0 )
        TGraphXclusters1->Draw("P");
    if ( TGraphXclusters2->GetN() > 0 )
        TGraphXclusters2->Draw("P");
    anXtrack = Reconstruct(TGraphXclusters2);
    anXtrack.Draw();
  
    myEventDisplay->cd(2);
    TGraphYhits->Draw("P");
    /*
      if ( TGraphYclusters->GetN() > 0 )
      TGraphYclusters->Draw("P");
      anYtrack = Reconstruct(TGraphYclusters);
    */
    if ( TGraphYclusters1->GetN() > 0 )
        TGraphYclusters1->Draw("P");
    if ( TGraphYclusters2->GetN() > 0 )
        TGraphYclusters2->Draw("P");
    anYtrack = Reconstruct(TGraphYclusters2);
    anYtrack.Draw();

    myEventDisplay->cd();
    myEventDisplay->Update();

    // print debug

    /*
      std::cout << "TGraphXhits: " << TGraphXhits->GetN() << std::endl;
      TGraphXhits->Print();
      std::cout << "TGraphYhits: " << TGraphYhits->GetN() << std::endl;
      TGraphYhits->Print();
      std::cout << "TGraphXclusters1: " << TGraphXclusters1->GetN() << std::endl;
      TGraphXclusters1->Print();
      std::cout << "TGraphYclusters1: " << TGraphYclusters1->GetN() << std::endl;
      TGraphYclusters1->Print();
      std::cout << "TGraphXclusters2: " << TGraphXclusters2->GetN() << std::endl;
      TGraphXclusters2->Print();
      std::cout << "TGraphYclusters2: " << TGraphYclusters2->GetN() << std::endl;
      TGraphYclusters2->Print();
      recon->PrintTkrTrk1Clusters();
      std::cout << "anXtrack" << std::endl;
      anXtrack.Print();
      std::cout << "anYtrack" << std::endl;
      anYtrack.Print();
    */
}

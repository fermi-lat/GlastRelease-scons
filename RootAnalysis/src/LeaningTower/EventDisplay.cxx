#include "EventDisplay.h"

#define TOWER 1

ClassImp(EventDisplay)

EventDisplay::EventDisplay(TString filename, TString geoFileName,int temid) 
  : m_temid(temid) {
    gStyle->SetCanvasColor(10);
    myTracker = new Tracker;
    if ( TOWER ) {
        if ( geoFileName.Length() == 0 )
            geoFileName = "$ROOTANALYSISROOT/src/LeaningTower/geometry/Tower0Geometry.txt";
        myTracker->SetTower(TOWER);
    }
    else
        geoFileName = "$ROOTANALYSISROOT/src/LeaningTower/geometry/stack2geometry.txt";
    myTracker->loadGeometry(geoFileName);

    myEvent = new Event(filename, myTracker->GetGeometry());
  
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

    TPaveStat = new TPaveText(300, 710, 500, 785);
    TPaveStat->SetTextAlign(12);
    TPaveStat->SetTextSize(0.025);

    runId = -1;
}

void EventDisplay::Go(int numEvent) {
    // numEvent is the number of the record in the tree, not the event id
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
    int running=0;
    do{
      myEvent->Go(entry+running);
      running++;
    } while(myEvent->GetTemId()!=m_temid);
    --running;
    entry+=running;
    std::cout<<"Passed "<<running<<" events"<<std::endl;
    myEventDisplay->cd();

    if ( runId == -1 )
        runId = myEvent->GetRunId();
    const int eventId         = myEvent->GetEventId();
    const float ebfTime       = myEvent->GetEbfTime();
    const int TkrTotalNumHits = myEvent->GetTkrTotalNumHits();

    std::cout << "RunId           = " << runId << std::endl;
    std::cout << "EventId         = " << eventId << std::endl;
    std::cout << "EbfTime         = " << ebfTime << std::endl;
    std::cout << "TkrTotalNumHits = " << TkrTotalNumHits << std::endl;

    TPaveStat->Clear();
    TPaveStat->AddText((TString("RunId:     "))+=runId);
    TPaveStat->AddText((TString("EventId:  "))+=eventId);
    TPaveStat->AddText((TString("EbfTime: "))+=ebfTime);
    myEventDisplay->cd(2);
    TPaveStat->Draw();

    Recon* recon = myEvent->GetRecon();

    TList *myGeometry = myTracker->GetGeometry();
    // this here "corrects" cluster positions with respect to the recon results.
    recon->TkrAlignmentSvc(myGeometry);

    int TkrNumClus = recon->GetTkrNumClus();
    int TkrNumTracks = recon->GetTkrNumTracks();
    int TkrTrk1NumClus = recon->GetTkrTrk1NumClus();
    std::cout << "TkrNumClus      = " << TkrNumClus << std::endl;
    std::cout << "TkrNumTracks    = " << TkrNumTracks << std::endl;
    std::cout << "TkrTrk1NumClus  = " << TkrTrk1NumClus << std::endl;
    if ( TkrTrk1NumClus > TkrNumClus || TkrNumClus > TkrTotalNumHits )
        std::cout << "??? TkrTrk1NumClus > TkrNumClus > TkrTotalNumHits ???"
                  << std::endl;

    TIter next(myGeometry);
    int i=0;
    while ( Layer* aPlane = (Layer*)next() ) {
        const double height = aPlane->GetHeight();
        const char* planeName = aPlane->GetName();
        const int planeNumHits = myEvent->GetPlaneNumHits(planeName);
        const int ToT0 = myEvent->GetToT(planeName, 0);
        const int ToT1 = myEvent->GetToT(planeName, 1);

        if ( aPlane->IsX() ) 
            myEventDisplay->cd(1);      
        else
            myEventDisplay->cd(2);

        TString dummy;
        if ( ToT0 >= 0 )
            dummy += ToT0;
        dummy += aPlane->GetTriggerReq(false) ? " x" : " .";
        TriggerReqText[i][0].SetText(-5, height, dummy);
        dummy = aPlane->GetTriggerReq(true) ? "x " : ". ";
        if ( ToT1 >= 0 )
            dummy += ToT1;
        TriggerReqText[i][1].SetText(365, height, dummy);
        LabelNumHits[i].SetText(425, height, (TString("(")+=planeNumHits)+")");
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
    myEventDisplay->SetTitle((TString("#")+=eventId));
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

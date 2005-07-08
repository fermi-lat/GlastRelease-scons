#include "Efficiency.h"

ClassImp(Efficiency)

Efficiency::Efficiency(const TString filename, const TString effFileName,
                       TString geoFileName, bool d) : debug(d) {
    if ( geoFileName.Length() == 0 )
        geoFileName = "$ROOTANALYSISROOT/src/LeaningTower/geometry/TowerBgeometry306000517.txt";
    myTracker = new Tracker;
    myTracker->loadGeometry(geoFileName);
    myTracker->SetTower(true);
    myEvent = new Event(filename, myTracker->GetGeometry());
    myEffFileName = effFileName;
}

void Efficiency::Go(int lastEntry) {
    int numEntries = myEvent->GetEntries();
    if ( lastEntry == -1 )
        lastEntry = numEntries;
    lastEntry = std::min(numEntries, lastEntry);

    TFile f(myEffFileName, "recreate");
    TTree effTree("efficiencyTree", "efficiencyTree");
    struct Ntuple ntuple0;
    effTree.Branch("eventId", &ntuple0.eventId, "event id/I");
    effTree.Branch("plane", ntuple0.charPlaneName, "name of plane[4]/C");
    effTree.Branch("xExt", &ntuple0.xExt, "extrapolated x position of hit/F");
    effTree.Branch("yExt", &ntuple0.yExt, "extrapolated y position of hit/F");
    effTree.Branch("siDist", &ntuple0.siDist, "distance to active area/F");
    effTree.Branch("res", &ntuple0.res,
                   "res (ext - data) to the closest cluster/F");
    effTree.Branch("siDistOther", &ntuple0.siDistOther,
                   "distance to active area in the other view/F");
    effTree.Branch("resOther", &ntuple0.resOther,
                 "res (ext - data) to the closest cluster in the other view/F");
    effTree.Branch("trigger", &ntuple0.trigger,
                 "event would have triggered even without this plane/B");

    const TList* myGeometry = myTracker->GetGeometry();
    const std::vector<TString> planeCols[2] = {
        myTracker->GetPlaneNameCol(0), myTracker->GetPlaneNameCol(1) };

    Recon* recon = myEvent->GetRecon();
    Progress progress;

    // some counters for the statistics
    int usedTrack[2] = { 0, 0 };
    int usedEvent = 0;

    for ( int entry=0; entry<lastEntry; ++entry ) { 
        myEvent->Go(entry);
        const Int_t eventId = myEvent->GetEventId();

        progress.Go(entry, lastEntry);

        recon->GetEvent(entry);
        //        const Int_t TkrNumClus = recon->GetTkrNumClus();
        //        const Int_t TkrTrk1NumClus = recon->GetTkrTrk1NumClus();

        // boy, this cut killed me.
        // With the new TkrRecon (or my omission of some jobOptions parameter),
        // the broken planes in tower A cause tracks to split into two.
        // Result: these events are not taken, and raise the efficiency.
        /*
        const Int_t TkrNumTracks = recon->GetTkrNumTracks();
        if ( TkrNumTracks != 1 )
            continue;
        */

        // "corrects" cluster positions with respect to the recon results
        recon->TkrAlignmentSvc(myGeometry);

        TLine TLtracks[2];
        for ( int view=0; view<2; ++view ) {
            // adding all clusters on planeCols[view] planes to TGclus
            const TGraph TGclusters = recon->GetClusterGraph(planeCols[view]);
            TLtracks[view] = Reconstruct(&TGclusters);
        }

        if( TLtracks[0].GetLineStyle() != 1 || TLtracks[1].GetLineStyle() != 1 )
            continue;   // no "good" Chi2

        // helper for the statistics
        //        bool trackUsed[2] = { false, false };
        //        bool eventUsed = false;

        // with this cut, the counters are redundant
        ++usedEvent;
        ++usedTrack[0];
        ++usedTrack[1];

        // I have to require that we have layers, i.e. pairs of xy planes
        // I use the event id as a flag to check that all planes got filled
        Ntuple ntuple[18][2];
        for ( int l=0; l<18; ++l )
            for ( int v=0; v<2; ++v )
                ntuple[l][v].eventId = -1;
        // now, extrapolate the fitted tracks in both views over all planes
        TIter next(myGeometry);
        while ( Layer* plane = (Layer*)next() ) {
            const TString planeName = plane->GetName();
            const int l = plane->GetLayer();
            const int v = plane->GetView();
            ntuple[l][v].eventId = eventId;
            strncpy(ntuple[l][v].charPlaneName, planeName.Data(), 4);
            const float posZ = plane->GetHeight();
            // x and y here are absolute, not in the system of the plane
            ntuple[l][v].xExt = ExtrapolateCoordinate(TLtracks[0], posZ);
            ntuple[l][v].yExt = ExtrapolateCoordinate(TLtracks[1], posZ);
            float absExt, ordExt;
            if ( plane->GetView() ) { // Y
                absExt = ntuple[l][v].yExt;
                ordExt = ntuple[l][v].xExt;
            }
            else { // X
                absExt = ntuple[l][v].xExt;
                ordExt = ntuple[l][v].yExt;
            }

            ntuple[l][v].siDist = plane->activeAreaDist(absExt, ordExt);

            ntuple[l][v].res = 400; // default if no clusters in this plane
            TGraph clusters = recon->GetClusterGraph(planeName);
            for ( int i=0; i<clusters.GetN(); ++i ) {
                Double_t hor, vert;
                clusters.GetPoint(i, hor, vert);
                const float newRes = absExt - hor;
                if ( std::abs(newRes) < std::abs(ntuple[l][v].res) )
                     ntuple[l][v].res = newRes;
            }
            // check if this plane is needed for a 3-in-a-row trigger
            if ( Get_3_in_a_row(l) )
                ntuple[l][v].trigger = true;
            else
                ntuple[l][v].trigger = false;
        }
        for ( int l=0; l<18; ++l )
            for ( int v=0; v<2; ++v ) {
                // check first that the ntuple was filled
                if ( ntuple[l][v].eventId < 0 )
                    std::cerr << "in layer " << l << " view " << v
                              << ": event id less than 0: "
                              << ntuple[l][v].eventId << std::endl;
                ntuple0 = ntuple[l][v];
                const int ov = v ? 0 : 1;
                ntuple0.siDistOther = ntuple[l][ov].siDist;
                ntuple0.resOther = ntuple[l][ov].res;
                effTree.Fill();
            }
    }

    std::cout << "used " << usedEvent << " of " << lastEntry << std::endl;
    std::cout << "found " << ++usedTrack[0] << " good tracks in X and "
              << ++usedTrack[1] << " in Y" << std::endl;

    f.Write();
    f.Close();
}

bool Efficiency::Get_3_in_a_row(const int inThisLayer) const {
    // let's first check if a single layer would have triggered
    bool layer[18];
    for ( int l=0; l<18; ++l ) {
	const TString nameX = GetPlaneName(l, 0);
	const TString nameY = GetPlaneName(l, 1);
	// OR of left and right for each view, and AND of both views
	layer[l] = ( myEvent->GetTriggerReq(nameX, 0)
                     || myEvent->GetTriggerReq(nameX, 1) )
	    && ( myEvent->GetTriggerReq(nameY, 0)
                 || myEvent->GetTriggerReq(nameY, 1) );
    }
    layer[inThisLayer] = false;  // in this layer is the plane we want to study!

    // do we still have a 3-in-a-row?
    for ( int l=0; l<16; ++l )
	if ( layer[l] && layer[l+1] && layer[l+2] )
	    return true;
    return false;
}
  
void Efficiency::Draw2D(const TString planeName, TCut cut,
                        const float residualDist, const float borderWidth)const{
    if ( planeName != "all" )
        cut += isPlane(planeName);

    TString canvasName = "Draw2D";
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c )
        c = new TCanvas(canvasName, canvasName, 800, 800);
    c->cd();
    gPad->SetTicks(1,1);
    //    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    TH1F* hframe = c->DrawFrame(-20, -20, 375, 375);
    hframe->GetXaxis()->SetTitle("x/mm");
    hframe->GetYaxis()->SetTitle("y/mm");
    hframe->SetTitle(planeName);
    c->Update();

    TFile f(myEffFileName, "READ");
    TTree* t = (TTree*)f.Get("efficiencyTree");
    t->SetMarkerSize(0.7);

    // draw all hits found in the active area:
    // ... not close to the edge
    t->SetMarkerColor(8);  // green
    t->Draw("yExt:xExt",
            cut && hitFoundDeepInActiveArea(residualDist,borderWidth), "same");
    c->Update();
    // ... close to the edge
    t->SetMarkerColor(6);  // magenta
    t->Draw("yExt:xExt", cut && hitFoundCloseToEdge(residualDist,borderWidth),
            "same");
    c->Update();

   // draw all interceptions which are not in an active area
    t->SetMarkerColor(2);  // red
    t->SetMarkerStyle(1);
    t->Draw("yExt:xExt", cut&&!inActiveArea(), "same");
    c->Update();

    // draw all missing hits in the active area:
    // ... not close to the edge
    t->SetMarkerColor(1);   // black
    t->SetMarkerStyle(28);  // swiss cross
    t->Draw("yExt:xExt",
            cut&&missingHitDeepInActiveArea(residualDist,borderWidth), "same");
    c->Update();
    // ... close to the edge
    t->SetMarkerColor(16);  // gray
    t->SetMarkerStyle(28);  // swiss cross
    t->Draw("yExt:xExt", cut && missingHitCloseToEdge(residualDist,borderWidth),
            "same");
    c->Update();

    // draw all interceptions with the active area which are neither hits nor
    // missing hits (i.e. rejected missing hits)
    t->SetMarkerColor(1);  // black
    t->SetMarkerStyle(5);  // x
    t->Draw("yExt:xExt", cut&&rejectedMissingHit(residualDist), "same");
    c->Update();

    if ( debug ) {
        t->SetScanField(0);
        t->Scan("eventId:plane:xExt:yExt",
                cut && inActiveArea()&&missingHit(residualDist));
        /*
        t->SetMarkerColor(2);  // red
        t->SetMarkerStyle(4);  // o
        t->SetMarkerSize(1.5);
        std::vector<int> eventCol;
        eventCol.push_back(343);
        ...
        TCut eventCut;
        for ( std::vector<int>::iterator it = eventCol.begin();
              it != eventCol.end(); ++it )
            eventCut = eventCut || (TString("eventId==")+=(*it));
        eventCut.Print();
        t->Draw("yExt:xExt", eventCut, "same");
        c->Update();
        */
    }

    f.Close();
}

void Efficiency::GetEfficiency(const TString planeName, const TCut cut,
                               const float residualDist,
                               const float borderWidth) const {

    TFile f(myEffFileName, "READ");
    TTree* t = (TTree*)f.Get("efficiencyTree");

    cut.Print();
    std::cout << " Object   Plane/L/W  Efficiency  Inefficiency     "
              << "hits  missed     all" << std::endl;

    const bool makeSummary = planeName == "all" || planeName == "plane";
    const TCut isMissingHit =
        missingHitDeepInActiveArea(residualDist, borderWidth);
    const TCut isHit(hitFoundDeepInActiveArea(residualDist, borderWidth));
    const TList* myGeometry = myTracker->GetGeometry();
    TIter next(myGeometry);
    while ( Layer* thePlane = (Layer*)next() ) {
        const TString thePlaneName = thePlane->GetName();
        if ( makeSummary || thePlaneName == planeName ) {
            const TCut planeCut(cut&&isPlane(thePlaneName));
            if ( planeName != "plane" ) {
                const TString ladderPar = thePlane->GetView() ? "yExt" : "xExt";
                const TString waferPar  = thePlane->GetView() ? "xExt" : "yExt";
                for ( int ladder=0; ladder<4; ++ladder ) {
                    TString ladderName(TString()+=ladder);
                    // the cuts should come from geometry in class Layer
                    // dirty hack, as always preliminary
                    TCut xmin((TString(ladderPar+">=")+=ladder*89.5).Data());
                    TCut xmax((
                               TString(ladderPar+"<=") += (ladder+1) * 89.5
                               ).Data());
                    const TCut ladderCut(planeCut&&xmin&&xmax);
                    for ( int wafer=0; wafer<4; ++wafer ) {
                        const TString waferName(TString()+=wafer);
                        TCut ymin((TString(waferPar+">=")+=wafer*89.5).Data());
                        TCut ymax((
                                   TString(waferPar+"<=") += (wafer+1) * 89.5
                                   ).Data());
                        const TCut waferCut(ladderCut&&ymin&&ymax);
                        const Long64_t numMiss =
                            t->Draw("eventId", waferCut&&isMissingHit, "goff");
                        PrintEfficiency(thePlaneName+'/'+ladderName+'/'
                                        +waferName,
                              t->Draw("eventId",waferCut&&isHit,"goff")+numMiss,
                                        numMiss);
                    }
                    const Long64_t numMiss =
                        t->Draw("eventId", ladderCut&&isMissingHit, "goff");
                    PrintEfficiency(thePlaneName+'/'+ladderName+"  ",
                             t->Draw("eventId",ladderCut&&isHit,"goff")+numMiss,
                                    numMiss);
                }
            }
            const Long64_t numMiss =
                t->Draw("eventId", planeCut&&isMissingHit, "goff");
            PrintEfficiency(thePlaneName+"    ",
                            t->Draw("eventId",planeCut&&isHit,"goff") + numMiss,
                            numMiss, t->Draw("eventId", planeCut, "goff"));
        }
    }
    if ( makeSummary ) {
        const Long64_t numMiss = t->Draw("eventId", cut&&isMissingHit, "goff");
        PrintEfficiency("           ",
                        t->Draw("eventId", cut&&isHit, "goff") + numMiss,
                        numMiss, t->Draw("eventId", cut, "goff"));
    }

    f.Close();
}

void Efficiency::PrintEfficiency(TString planeName,  const int hits,
                                 const int missing, const int all)
    const {
    for ( ; planeName.Length()<11; )
                planeName.Prepend(' '); 
    std::cout << std::setw(7);
    if ( planeName(10) != ' ' )
        std::cout << "  wafer";
    else if ( planeName(8) != ' ' )
        std::cout << " ladder";
    else if ( planeName(6) != ' ' )
        std::cout << "  plane";
    else
        std::cout << "  tower";
    const float ineff = hits ? 100. * missing / hits : -100;
    std::cout << ' ' << std::setw(11) << planeName << ' '
              << std::setiosflags(std::ios::fixed)
              << std::setw(9) << std::setprecision(2) << 100. - ineff << " % "
              << std::setw(11) << std::setprecision(2) << ineff << " % "
              << std::setw(8) << hits
              << std::setw(8) << missing;
    if ( all )
        std::cout << std::setw(8) << all;
    std::cout << std::endl;
}

void Efficiency::DrawEfficiency(const TString planeName, TCut cut,
                                const float residualDist,
                                const float borderWidth) const{
    TFile f(myEffFileName, "READ");
    TTree* t = (TTree*)f.Get("efficiencyTree");
    gROOT->cd();

    cut += isPlane(planeName);

    Layer* thePlane = (Layer*)myTracker->GetGeometry()->FindObject(planeName);
    if ( !thePlane ) {
        std::cerr << "plane " << planeName << " not found!" << std::endl;
        return;
    }

    const Int_t nBins(100);
    const Double_t xmin = -20.;
    const Double_t xmax = 380.;
    const TString allhitsName("allHits"+planeName);
    TH1F* hAll = (TH1F*)gROOT->FindObject(allhitsName);
    if ( !hAll ) {
        hAll = new TH1F(allhitsName, allhitsName, nBins, xmin, xmax);
        hAll->Sumw2();
    }
    hAll->GetXaxis()->SetTitle("pos/mm");
    hAll->GetYaxis()->SetTitle("num");
    const TString var(thePlane->GetView()?"yExt":"xExt");
    t->Project(allhitsName, var,
         cut&&deepInActiveArea(borderWidth)&&!rejectedMissingHit(residualDist));

    const TString missName("missingHits"+planeName);
    TH1F* hBad = (TH1F*)gROOT->FindObject(missName);
    if ( !hBad ) {
        hBad = new TH1F(missName, missName, nBins, xmin, xmax);
        hBad->Sumw2();
    }
    hBad->GetXaxis()->SetTitle("pos/mm");
    hBad->GetYaxis()->SetTitle("num");
    t->Project(missName, var,
            cut&&missingHitDeepInActiveArea(residualDist, borderWidth), "hist");

    const TString ineffName("inefficiency"+planeName);
    TH1F* hIneff = (TH1F*)gROOT->FindObject(ineffName);
    if ( !hIneff )
        hIneff = new TH1F(ineffName, ineffName, nBins, xmin, xmax);
    hIneff->GetXaxis()->SetTitle("pos/mm");
    hIneff->GetYaxis()->SetTitle("inefficiency");
    hIneff->Divide(hBad, hAll);

    f.Close();

    // here, do the drawing
    const TString canvasName = "DrawEfficiency" + planeName;
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c ) {
        c = new TCanvas(canvasName, canvasName, 600, 800);
        c->Divide(1,3);
    }
    c->cd();

    c->cd(1);
    gPad->SetTicks(1,1);
    hAll->Draw("hist");

    c->cd(2);
    gPad->SetTicks(1,1);
    hBad->Draw("hist");

    c->cd(3);
    gPad->SetTicks(1,1);
    hIneff->Draw("e0");

    c->cd();
}

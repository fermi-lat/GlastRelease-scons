#include "TCut.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLine.h"
#include "TList.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "Tracker.h"
#include "Layer.h"
#include "Event.h"
#include "Recon.h"
#include "Progress.h"

class Efficiency {
private:
    Event*   myEvent;
    Tracker* myTracker;
    TString  myEffFileName;
public:
    Efficiency(const TString="MyRootFile.root", const TString="efficiency.root",
               TString geoFileName="");
    ~Efficiency() {
        delete myEvent;
        delete myTracker;
    }
    // structure to save the ntuple values before filling the tree
    struct Ntuple {
        Int_t eventId;
        Char_t charPlaneName[4];
        Float_t xExt, yExt, siDist, res;
        Float_t siDistOther, resOther;
    };

    void Go(int lastEntry=-1);

    void Draw2D(TString plane="all", TCut="", float residualDist=1,
                float borderWidth=1) const;
    void GetEfficiency(const TString plane="all", const TCut="", const float residualDist=1,
                       const float borderWidth=1, const bool draw=false) const;
    void DrawEfficiency(TString plane, TCut="", float residualDist=1,
                        float borderWidth=1, bool draw=true)
        const;
private:
    void PrintEfficiency(TString plane, float ineff, int hits, int missing) const;
};

Efficiency::Efficiency(const TString filename, const TString effFileName,
                       TString geoFileName) {
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
    effTree.Branch("res", &ntuple0.res, "res (ext - data) to the closest cluster/F");
    effTree.Branch("siDistOther", &ntuple0.siDistOther, "distance to active area in the other view/F");
    effTree.Branch("resOther", &ntuple0.resOther, "res (ext - data) to the closest cluster in the other view/F");

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
        const Int_t TkrNumTracks = recon->GetTkrNumTracks();
        //        const Int_t TkrTrk1NumClus = recon->GetTkrTrk1NumClus();

        if ( TkrNumTracks != 1 )
            continue;

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
        // I use the event id as a flag
        Ntuple ntuple[18][2];
        for ( int l=0; l<18; ++l )
            for ( int v=0; v<2; ++v )
                ntuple[l][v].eventId = -1;
        // now, loop over all planes
        TIter next(myGeometry);
        while ( Layer* plane = (Layer*)next() ) {
            const TString planeName = plane->GetName();
            const int ll = plane->GetLayer();
            const int vv = plane->GetView();
            ntuple[ll][vv].eventId = eventId;
            strncpy(ntuple[ll][vv].charPlaneName, planeName.Data(), 4);
            const float posZ = plane->GetHeight();
            // x and y here are absolute, not in the system of the plane
            ntuple[ll][vv].xExt = ExtrapolateCoordinate(TLtracks[0], posZ);
            ntuple[ll][vv].yExt = ExtrapolateCoordinate(TLtracks[1], posZ);
            float absExt, ordExt;
            if ( plane->GetView() ) { // Y
                absExt = ntuple[ll][vv].yExt;
                ordExt = ntuple[ll][vv].xExt;
            }
            else { // X
                absExt = ntuple[ll][vv].xExt;
                ordExt = ntuple[ll][vv].yExt;
            }

            ntuple[ll][vv].siDist = plane->activeAreaDist(absExt, ordExt);

            ntuple[ll][vv].res = 400;
            TGraph clusters = recon->GetClusterGraph(planeName);
            for ( int i=0; i<clusters.GetN(); ++i ) {
                Double_t h, v;
                clusters.GetPoint(i, h, v);
                const float newRes = absExt - h;
                if ( std::abs(newRes) < std::abs(ntuple[ll][vv].res) )
                     ntuple[ll][vv].res = newRes;
            }
            //            effTree.Fill();
        }
        // it shouldn't happen, so lets check first if all ntuple were filled
        for ( int l=0; l<18; ++l )
            for ( int v=0; v<2; ++v )
                if ( ntuple[l][v].eventId < 0 )
                    std::cerr << "in layer " << l << " view " << v
                              << ": event id less than 0: "
                              << ntuple[l][v].eventId << std::endl;
        for ( int l=0; l<18; ++l )
            for ( int v=0; v<2; ++v ) {
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

void Efficiency::Draw2D(TString plane, TCut cut, float residualDist,
                        float borderWidth) const {
    if ( plane != "all" ) {
        const TCut cutPlane("plane==\"" + plane + "\"");
        cut += cutPlane;
    }

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
    hframe->SetTitle(plane);
    c->Update();

    TFile f(myEffFileName, "READ");
    TTree* t = (TTree*)f.Get("efficiencyTree");
    t->SetMarkerSize(0.7);

    TCut inActiveArea("siDist<0");
    TCut closeToEdge(inActiveArea+(TString("siDist>-")+=borderWidth).Data());
    TCut hitFound((TString("res<")+=residualDist).Data());
    // cuts in the other view
    //    TCut inActiveAreaOther("siDistOther<0");
    //    TCut closeToEdgeOther(inActiveAreaOther+(TString("siDistOther>-")+=borderWidth).Data());
    TCut hitFoundOther((TString("resOther<")+=residualDist).Data());
    // missingHit is != !hitFound
    // a hit is missing, if:
    //   1) we don't find a hit in the plane within hitFound
    //   2) we find a hit in the other plane of the same layer within hitFoundOther.
    // If also in the other plane the hit is missing, it is very probable that the particle got stopped,
    // but we continue tracing the track through the tower.  Hopefully, we don't get screwed by noise hits.
    TCut missingHit = !hitFound && hitFoundOther;

   // draw all interceptions which are not in an active area
    t->SetMarkerColor(2);
    t->SetMarkerStyle(1);
    t->Draw("yExt:xExt", cut&&!inActiveArea, "same");
    c->Update();

    // draw all hits in the active area close to the edge
    t->SetMarkerColor(6);
    t->Draw("yExt:xExt", cut&&closeToEdge&&hitFound, "same");
    c->Update();

    // draw all hits in the active area not close to the edge
    t->SetMarkerColor(8);
    t->Draw("yExt:xExt", cut&&inActiveArea&&!closeToEdge&&hitFound, "same");
    c->Update();

    // draw all missing hits in the active area but close to the edge
    t->SetMarkerColor(16);
    t->SetMarkerStyle(28);
    t->Draw("yExt:xExt", cut&&closeToEdge&&missingHit, "same");
    c->Update();

    // draw all missing hits in the active area not close to the edge
    t->SetMarkerColor(1);
    t->SetMarkerStyle(28);
    t->Draw("yExt:xExt", cut&&inActiveArea&&!closeToEdge&&missingHit, "same");
    c->Update();

    // draw all interceptions with the active area which are not hits and not missing hits
    t->SetMarkerColor(1);
    t->SetMarkerStyle(5);
    t->Draw("yExt:xExt", cut&&inActiveArea&&!hitFound&&!missingHit, "same");
    c->Update();

    f.Close();
}

void Efficiency::GetEfficiency(const TString planeName, const TCut cut, const float residualDist,
                               const float borderWidth, const bool draw) const {
    cut.Print();
    std::cout << " Object   Plane Ladder  Wafer  Efficiency  Inefficiency     Hits  Missing hits" << std::endl;

    const bool makeSummary = planeName == "all" || planeName == "plane";
    int totalHits = 0;
    int totalMissing = 0;
    const TList* myGeometry = myTracker->GetGeometry();
    TIter next(myGeometry);
    while ( Layer* thePlane = (Layer*)next() ) {
        const TString thePlaneName = thePlane->GetName();
        if ( makeSummary || thePlaneName == planeName ) {
            TString name(thePlaneName);
            for ( ; name.Length()<7; )
                name.Prepend(' ');
            if ( planeName != "plane" ) {
                TString ladderPar = thePlane->GetView() ? "yExt" : "xExt";
                TString waferPar  = thePlane->GetView() ? "xExt" : "yExt";
                for ( int ladder=0; ladder<4; ++ladder ) {
                    TString ladderName(TString()+=ladder);
                    for ( ; ladderName.Length()<7; )
                        ladderName.Prepend(' ');
                    // the cuts should come from geometry in class Layer
                    // dirty hack, as always preliminary

                    TCut xmin((TString(ladderPar+">=")+=ladder*89.5).Data());
                    TCut xmax((TString(ladderPar+"<=")+=(ladder+1)*89.5).Data());
                    TCut ladderCut = cut && xmin && xmax;
                    for ( int wafer=0; wafer<4; ++wafer ) {
                        TString waferName(TString()+=wafer);
                        for ( ; waferName.Length()<7; )
                            waferName.Prepend(' ');
                        TCut ymin((TString(waferPar+">=")+=wafer*89.5).Data());
                        TCut ymax((TString(waferPar+"<=")+=(wafer+1)*89.5).Data());
                        TCut waferCut = ladderCut && ymin && ymax;
                        DrawEfficiency(thePlaneName, waferCut, residualDist, borderWidth, draw);
                        PrintEfficiency(name+ladderName+waferName, 100.*thePlane->GetInefficiency(),
                                        thePlane->GetHitsInActiveArea(), thePlane->GetMissingHits());
                    }
                    DrawEfficiency(thePlaneName, ladderCut, residualDist, borderWidth, draw);
                    PrintEfficiency(name+ladderName+"       ", 100.*thePlane->GetInefficiency(),
                                    thePlane->GetHitsInActiveArea(), thePlane->GetMissingHits());
                }
            }
            DrawEfficiency(thePlaneName, cut, residualDist, borderWidth, draw);
            PrintEfficiency(name+"              ", 100.*thePlane->GetInefficiency(),
                            thePlane->GetHitsInActiveArea(), thePlane->GetMissingHits());
            if ( makeSummary ) {
                totalHits += thePlane->GetHitsInActiveArea();
                totalMissing += thePlane->GetMissingHits();
            }
        }
    }
    if ( makeSummary )
        PrintEfficiency("                     ", 100.*totalMissing/totalHits, totalHits, totalMissing);
}

void Efficiency::PrintEfficiency(TString planeName, float ineff, int hits, int missing)
    const {
    TString object;
    if ( planeName(20) != ' ' )
        object = "  wafer";
    else if ( planeName(13) != ' ' )
        object = " ladder";
    else if ( planeName(6) != ' ' )
        object = "  plane";
    else
        object = "  tower";
    std::cout << std::setw(7) << object << ' '
              << std::setw(21) << planeName << ' '
              << std::setiosflags(std::ios::fixed)
              << std::setw(9) << std::setprecision(1) << 100. - ineff << " % "
              << std::setw(11) << std::setprecision(1) << ineff << " % "
              << std::setw(8) << hits
              << std::setw(14) << missing
              << std::endl;
}

void Efficiency::DrawEfficiency(TString planeName, TCut cut, float residualDist,
                                float borderWidth, bool draw) const {
    TFile f(myEffFileName, "READ");
    TTree* t = (TTree*)f.Get("efficiencyTree");
    gROOT->cd();

    const TCut cutPlane("plane==\"" + planeName + "\"");
    cut += cutPlane;

    TString var;
    Layer* thePlane = (Layer*)myTracker->GetGeometry()->FindObject(planeName);
    if ( !thePlane ) {
        std::cerr << "plane " << planeName << " not found!" << std::endl;
        return;
    }
    if ( thePlane->GetView() )
        var = "yExt";
    else
        var = "xExt";

    TCut inActiveArea("siDist<0");
    TCut closeToEdge(inActiveArea+(TString("siDist>-")+=borderWidth).Data());
    TCut hitFound((TString("res<")+=residualDist).Data());
    TCut hitFoundOther((TString("resOther<")+=residualDist).Data());
    TCut missingHit(!hitFound&&hitFoundOther);
    TCut allHits = hitFound || missingHit;

    const Int_t nBins = 100;
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
    t->Project(allhitsName, var, cut&&inActiveArea&&!closeToEdge&&allHits);

    const TString missName("missingHits"+planeName);
    TH1F* hBad = (TH1F*)gROOT->FindObject(missName);
    if ( !hBad ) {
        hBad = new TH1F(missName, missName, nBins, xmin, xmax);
        hBad->Sumw2();
    }
    hBad->GetXaxis()->SetTitle("pos/mm");
    hBad->GetYaxis()->SetTitle("num");
    t->Project(missName, var, cut&&inActiveArea&&!closeToEdge&&missingHit, "hist");

    const TString ineffName("inefficiency"+planeName);
    TH1F* hIneff = (TH1F*)gROOT->FindObject(ineffName);
    if ( !hIneff )
        hIneff = new TH1F(ineffName, ineffName, nBins, xmin, xmax);
    hIneff->GetXaxis()->SetTitle("pos/mm");
    hIneff->GetYaxis()->SetTitle("inefficiency");
    hIneff->Divide(hBad, hAll);

    f.Close();

    thePlane->SetHitsInActiveArea(static_cast<int>(hAll->GetEntries()));
    thePlane->SetMissingHits(static_cast<int>(hBad->GetEntries()));

    if ( !draw )
        return;
    // here, doing the drawing
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

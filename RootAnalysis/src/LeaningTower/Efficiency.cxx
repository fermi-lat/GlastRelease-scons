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
    Efficiency(const TString="MyRootFile.root",const TString="efficiency.root");
    ~Efficiency() {
        delete myEvent;
        delete myTracker;
    }

    void Go(int lastEntry=-1);

    void Draw2D(TString plane="All", TCut="", float residualSize=0.228) const;
    void GetEfficiency(TString plane="All", TCut="", float residualSize=0.228,
                       bool draw=false) const;
    void PrintEfficiency(TString plane, float ineff, int missing) const;
    void DrawEfficiency(TString plane, TCut="", float residualSize=0.228,
                        bool draw=true) const;
};

Efficiency::Efficiency(TString filename, TString effFileName) {
    myTracker = new Tracker;
    myTracker->loadGeometry(gSystem->ExpandPathName(
             "$ROOTANALYSISROOT/src/LeaningTower/geometry/TowerAgeometry.txt"));
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
    Int_t eventId;
    Char_t charPlaneName[4];
    Float_t xExt, yExt, siDist, res;
    effTree.Branch("eventId", &eventId, "event id/I");
    effTree.Branch("plane", charPlaneName, "name of plane[4]/C");
    effTree.Branch("xExt", &xExt, "extrapolated x position of hit/F");
    effTree.Branch("yExt", &yExt, "extrapolated y position of hit/F");
    effTree.Branch("siDist", &siDist, "distance to active area/F");
    effTree.Branch("res", &res, "res (ext - data) to the closest cluster/F");

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
        eventId = myEvent->GetEventId();

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

        // now, loop over all planes
        TIter next(myGeometry);
        while ( Layer* plane = (Layer*)next() ) {
            const TString planeName = plane->GetName();

            strncpy(charPlaneName, planeName.Data(), 4);
            const float posZ = plane->GetHeight();
            // x and y here are absolute, not in the system of the plane
            xExt = ExtrapolateCoordinate(TLtracks[0], posZ);
            yExt = ExtrapolateCoordinate(TLtracks[1], posZ);
            float absExt, ordExt;
            if ( plane->GetView() ) { // Y
                absExt = yExt;
                ordExt = xExt;
            }
            else { // X
                absExt = xExt;
                ordExt = yExt;
            }

            siDist = plane->activeAreaDist(absExt, ordExt);

            res = 400;
            TGraph clusters = recon->GetClusterGraph(planeName);
            for ( int i=0; i<clusters.GetN(); ++i ) {
                Double_t h, v;
                clusters.GetPoint(i, h, v);
                const float newRes = absExt - h;
                if ( std::abs(newRes) < std::abs(res) )
                     res = newRes;
            }
            effTree.Fill();
        }
    }

    std::cout << "used " << usedEvent << " of " << lastEntry << std::endl;
    std::cout << "found " << ++usedTrack[0] << " good tracks in X and "
              << ++usedTrack[1] << " in Y" << std::endl;

    f.Write();
    f.Close();
}

void Efficiency::Draw2D(TString plane, TCut cut, float residualSize) const {
    if ( plane != "All" ) {
        const TCut cutPlane("plane==\"" + plane + "\"");
        cut += cutPlane;
    }

    TString canvasName = "Draw2D";
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c )
        c = new TCanvas(canvasName, canvasName, 800, 800);
    gPad->SetTicks(1,1);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    c->DrawFrame(-20, -20, 375, 375);

    TFile f(myEffFileName);
    TTree* t = (TTree*)f.Get("efficiencyTree");
    t->SetMarkerStyle(1);
    t->SetMarkerColor(1);

    TCut inActiveArea("siDist<0");
    TCut residualCut((TString("res>")+=residualSize).Data());

   // draw all interceptions which are not in an active area
    t->SetMarkerColor(2);
    t->Draw("yExt:xExt", cut+!inActiveArea, "same");

    // draw all interceptions with the active area (should be hits)
    t->SetMarkerColor(8);
    t->Draw("yExt:xExt", cut+inActiveArea, "same");

    // draw all missing hits in the active area
    t->SetMarkerColor(1);
    t->SetMarkerStyle(28);
    t->Draw("yExt:xExt", cut+inActiveArea+residualCut, "same");

    f.Close();
}

void Efficiency::GetEfficiency(TString planeName, TCut cut, float residualSize,
                               bool draw) const {
    cut.Print();
    std::cout << "residual>" << residualSize << std::endl;
    std::cout << "Plane  Efficiency  Inefficiency  Missing hits" << std::endl;

    int totalHits = 0;
    int totalMissing = 0;
    const TList* myGeometry = myTracker->GetGeometry();
    TIter next(myGeometry);
    while ( Layer* thePlane = (Layer*)next() ) {
        const TString thePlaneName = thePlane->GetName();
        if ( planeName == "All" || thePlaneName == planeName ) {
            DrawEfficiency(thePlaneName, cut, residualSize, draw);
            PrintEfficiency(thePlaneName, 100.*thePlane->GetInefficiency(),
                            thePlane->GetMissingHits());
            totalHits += thePlane->GetHitsInActiveArea();
            totalMissing += thePlane->GetMissingHits();
        }
    }
    PrintEfficiency("total", 100.*totalMissing/totalHits, totalMissing);
}

void Efficiency::PrintEfficiency(TString planeName, float ineff, int missing)
    const {
    std::cout << std::setw(5) << planeName << ' '
              << std::setiosflags(std::ios::fixed)
              << std::setw(9) << std::setprecision(1) << 100. - ineff << " % "
              << std::setw(11) << std::setprecision(1) << ineff << " % "
              << std::setw(13) << missing
              << std::endl;
}

void Efficiency::DrawEfficiency(TString planeName, TCut cut, float residualSize,
                                bool draw) const {
    TFile f(myEffFileName);
    TTree* t = (TTree*)f.Get("efficiencyTree");
    gROOT->cd();

    const TCut cutPlane("plane==\"" + planeName + "\"");
    const TCut inActiveArea("siDist<0");
    const TCut residualCut((TString("res>")+=residualSize).Data());
    cut += cutPlane + inActiveArea;

    TString var;
    Layer* thePlane = (Layer*)myTracker->GetGeometry()->FindObject(planeName);
    if ( thePlane->GetView() )
        var = "yExt";
    else
        var = "xExt";

    const Int_t nBins = 100;
    const Double_t xmin = -20.;
    const Double_t xmax = 380.;
    const TString goodName("activeHits"+planeName);
    TH1F* hGood = (TH1F*)gROOT->FindObject(goodName);
    if ( !hGood ) {
        hGood = new TH1F(goodName, goodName, nBins, xmin, xmax);
        hGood->Sumw2();
    }
    t->Project(goodName, var, cut);

    const TString badName("missingHits"+planeName);
    TH1F* hBad = (TH1F*)gROOT->FindObject(badName);
    if ( !hBad ) {
        hBad = new TH1F(badName, badName, nBins, xmin, xmax);
        hBad->Sumw2();
    }
    t->Project(badName, var, cut+residualCut, "hist");

    const TString ineffName("inefficiency"+planeName);
    TH1F* hIneff = (TH1F*)gROOT->FindObject(ineffName);
    if ( !hIneff )
        hIneff = new TH1F(ineffName, ineffName, nBins, xmin, xmax);
    hIneff->Divide(hBad, hGood);

    f.Close();

    thePlane->SetHitsInActiveArea(static_cast<int>(hGood->GetEntries()));
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

    c->cd(1);
    gPad->SetTicks(1,1);
    hGood->Draw("hist");

    c->cd(2);
    gPad->SetTicks(1,1);
    hBad->Draw("hist");

    c->cd(3);
    gPad->SetTicks(1,1);
    hIneff->Draw();

    c->cd();
}

#include "TFile.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMap.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "Tracker.h"
#include "Layer.h"
#include "Event.h"
#include "Recon.h"
#include "Progress.h"

class Residual {
private:
    Event*   myEvent;
    Tracker* myTracker;
public:
    Residual(const TString="MyRootFile.root");
    ~Residual() {
        delete myEvent;
        delete myTracker;
    }

    void Go(int lastEntry=-1);
};

Residual::Residual(TString filename) {
    myTracker = new Tracker;
    myTracker->loadGeometry(gSystem->ExpandPathName("$ROOTANALYSISROOT/src/LeaningTower/geometry/Tower0Geometry.txt"));
    myTracker->loadFitting(gSystem->ExpandPathName("$ROOTANALYSISROOT/src/LeaningTower/geometry/FittingPlanes.txt"));
    myTracker->IsTower(true);
    myEvent = new Event(filename, (TMap*)myTracker->GetGeometry());
}

void Residual::Go(int lastEntry) {
    int numEntries = myEvent->GetEntries();
    if ( lastEntry == -1 )
        lastEntry = numEntries;
    lastEntry = std::min(numEntries, lastEntry);

    int bins = 80000;
    double xmin = -400;
    double xmax = 400;
    TFile f("res.root", "recreate");

    TTree resTree("residualTree", "residualTree");
    Char_t charPlaneName[4];
    Float_t posHor, posVert, posHorExt, invSlope;
    resTree.Branch("name",     charPlaneName, "name of plane[4]/C");
    resTree.Branch("h",        &posHor,       "horizontal position of cluster/F");
    resTree.Branch("v",        &posVert,      "vertical position of cluster/F");
    resTree.Branch("hExt",     &posHorExt,    "extrapolated horizontal position of cluster/F");
    resTree.Branch("invSlope", &invSlope,     "inverse slope of track/F");

    TH1D Xresiduals("Xresiduals", "Xresiduals", bins, xmin, xmax);
    TH1D Yresiduals("Yresiduals", "Yresiduals", bins, xmin, xmax);

    TH1D hX[18];
    TH1D hY[18];
    for ( int i=0; i<18; ++i ) {
        TString title;
        title = TString("X") + (TString("residuals")+=i);
        hX[i] = TH1D(title, title, bins, xmin, xmax);
        title = TString("Y") + (TString("residuals")+=i);
        hY[i] = TH1D(title, title, bins, xmin, xmax);
    }

    //#define SANDWICH
#define ALLLAYERS
#ifdef SANDWICH
    const TMap *myGeometry = myTracker->GetGeometry();
#endif

    Recon* recon = myEvent->GetRecon();
    Progress progress;

    for ( int entry=0; entry<lastEntry; ++entry ) { 
        myEvent->Go(entry);

        progress.Go(entry, lastEntry);

        recon->GetEvent(entry);
        const Int_t TkrNumClus = recon->GetTkrNumClus();
        const Int_t TkrNumTracks = recon->GetTkrNumTracks();
        //        const Int_t TkrTrk1NumClus = recon->GetTkrTrk1NumClus();

        if ( TkrNumTracks != 1 )
            continue;

        // ATTENTION -------------------------------------------------
        // recon counts layers from top, i.e. X0 -> layer 17 view 0
        const Int_t* TkrClusLayer = recon->GetTkrClusLayer();
        const Int_t* TkrClusView = recon->GetTkrClusView();
        const Float_t* TkrClusX = recon->GetTkrClusX();
        const Float_t* TkrClusY = recon->GetTkrClusY();
        const Float_t* TkrClusZ = recon->GetTkrClusZ();

        for ( int i=0; i<TkrNumClus; ++i ) { // loop over all clusters
            const int iView = TkrClusView[i];
            const int iReconLayer = TkrClusLayer[i];
            const int iLayer = 17 - iReconLayer;
            const TString planeName = GetPlaneName(iLayer, iView);
#ifdef SANDWICH
            // version with sandwiching
            const Layer* plane = (Layer*)myGeometry->GetValue(planeName);
            const std::vector<TString> planeCol = plane->GetPlanesForFittingCol();
#endif
#ifdef ALLLAYERS
            // version with all layers
            std::vector<TString> planeCol = myTracker->GetPlaneNameCol(iView, true);  // make a list of "good" planes only!
            // remove all hits on the tray to which the plane to be studied belongs
            planeCol.erase(find(planeCol.begin(), planeCol.end(), planeName));
            planeCol.erase(find(planeCol.begin(), planeCol.end(), GetTwinPlaneName(iLayer, iView)));
#endif

            // adding all clusters on the planeCol layers to TGclus
            bool exactlyOne;
            const TGraph TGclus = recon->GetClusterGraph(planeCol, &exactlyOne);

            if ( !planeCol.size() || !exactlyOne )
                continue;
            // rejects if a bad plane (empty refCol), and if not each plane has exactly one cluster

            const TLine track = Reconstruct(&TGclus, false);

            if ( !IsValid(track) )
                continue;

            strncpy(charPlaneName, planeName.Data(), 4);
            posVert = TkrClusZ[i];
            posHorExt = ExtrapolateCoordinate(track, posVert);
            invSlope = ( track.GetX2() - track.GetX1() ) / ( track.GetY2() - track.GetY1() );
            if ( iView == 0 ) {
                posHor = TkrClusX[i];
                const double res = posHorExt - posHor;
                hX[iLayer].Fill(res);
                Xresiduals.Fill(res);
            }
            else if ( iView == 1 ) {
                posHor = TkrClusY[i];
                const double res = posHorExt - posHor;
                hY[iLayer].Fill(res);
                Yresiduals.Fill(res);
            }
            resTree.Fill();
        } // end of loop over all clusters
    } // end of loop over all events

    f.Write();
    f.Close();
}

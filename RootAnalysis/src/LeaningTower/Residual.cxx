#include "TCut.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMap.h"
#include "TString.h"
#include "TStyle.h"
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
    TString  myResFileName;
public:
    Residual(const TString="MyRootFile.root", const TString="residual.root");
    ~Residual() {
        delete myEvent;
        delete myTracker;
    }

    void Go(int lastEntry=-1);
    void DrawResidual(TCut="");
    void DrawResSlope(TCut="");
    void DrawResSlopeAll(TCut="abs(h_abs_ext-h_abs)<1");
    void DrawResOrd(TCut="");
    void DrawResOrdAll(TCut="abs(h_abs_ext-h_abs)<1");
};

Residual::Residual(TString filename, TString resFileName) {
    myTracker = new Tracker;
    myTracker->loadGeometry(gSystem->ExpandPathName(
             "$ROOTANALYSISROOT/src/LeaningTower/geometry/Tower0Geometry.txt"));
    myTracker->loadFitting(gSystem->ExpandPathName(
        "$ROOTANALYSISROOT/src/LeaningTower/geometry/Tower0FittingPlanes.txt"));
    myTracker->IsTower(true);
    myEvent = new Event(filename, (TMap*)myTracker->GetGeometry());
    myResFileName = resFileName;
}

void Residual::Go(int lastEntry) {
    int numEntries = myEvent->GetEntries();
    if ( lastEntry == -1 )
        lastEntry = numEntries;
    lastEntry = std::min(numEntries, lastEntry);

    int bins = 80000;
    double xmin = -400;
    double xmax = 400;
    TFile f(myResFileName, "recreate");

    TTree resTree("residualTree", "residualTree");
    Char_t charPlaneName[4];
    Float_t posHorAbs, posHorOrdExt, posVert, posHorAbsExt, invSlope;
    resTree.Branch("name", charPlaneName, "name of plane[4]/C");
    resTree.Branch("h_abs", &posHorAbs, "horizontal position of cluster/F");
    resTree.Branch("h_ord", &posHorOrdExt,
                   "from the other view: horizontal position of cluster/F");
    resTree.Branch("v", &posVert, "vertical position of cluster/F");
    resTree.Branch("h_abs_ext", &posHorAbsExt,
                   "extrapolated horizontal of cluster/F");
    resTree.Branch("invSlope", &invSlope, "inverse slope of track/F");

    TH1D hResiduals[2];
    hResiduals[0] = TH1D("Xresiduals", "Xresiduals", bins, xmin, xmax);
    hResiduals[1] = TH1D("Yresiduals", "Yresiduals", bins, xmin, xmax);

    TH1D hRes[2][18];
    for ( int i=0; i<18; ++i ) {
        TString title;
        title = TString("X") + (TString("residuals")+=i);
        hRes[0][i] = TH1D(title, title, bins, xmin, xmax);
        title = TString("Y") + (TString("residuals")+=i);
        hRes[1][i] = TH1D(title, title, bins, xmin, xmax);
    }

    const TMap* myGeometry = myTracker->GetGeometry();

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

        // "corrects" cluster positions with respect to the recon results
        recon->TkrAlignmentSvc(myGeometry);

        // ATTENTION -------------------------------------------------
        // recon counts layers from top, i.e. X0 -> layer 17 view 0
        const Int_t* TkrClusLayer = recon->GetTkrClusLayer();
        const Int_t* TkrClusView = recon->GetTkrClusView();
        const Float_t* TkrClusX = recon->GetTkrClusX();
        const Float_t* TkrClusY = recon->GetTkrClusY();
        const Float_t* TkrClusZ = recon->GetTkrClusZ();

        // lets first make two tracks which are later used to determine an
        // approximate position in the plane of interest.
        // We don't care that there is exactly one cluster per plane.  If it
        // would matter, the Chi2 would be bad anyway!
        TLine refTrack[2];
        for ( int i=0; i<2; ++i ) {
            // make a list of "good" planes only!
            const std::vector<TString> planeCol =
                myTracker->GetPlaneNameCol(i, true);
            // adding all clusters on the planeCol planes to TGclus
            const TGraph TGclus = recon->GetClusterGraph(planeCol);
            refTrack[i] = Reconstruct(&TGclus, false);
        }

        // now, loop over all clusters
        for ( int i=0; i<TkrNumClus; ++i ) {
            const int iView = TkrClusView[i];
            const int iReconLayer = TkrClusLayer[i];
            const int iLayer = 17 - iReconLayer;
            const TString planeName = GetPlaneName(iLayer, iView);
            std::vector<TString>planeCol=myTracker->GetPlaneNameCol(iView,true);
            // remove all hits on tray to which the plane to be studied belongs
            planeCol.erase(find(planeCol.begin(), planeCol.end(), planeName));
            planeCol.erase(find(planeCol.begin(), planeCol.end(),
                                GetTrayTwinPlaneName(iLayer, iView)));
            bool exactlyOne;
            const TGraph TGclus = recon->GetClusterGraph(planeCol, &exactlyOne);

            if ( !planeCol.size() || !exactlyOne )
                continue;
            // rejects if a bad plane (empty refCol), and if not each plane has
            // exactly one cluster

            const TLine track = Reconstruct(&TGclus, false);

            //            if ( !IsValid(track) ) // no track fitted
            if ( track.GetLineStyle() != 1 ) // no "good" Chi2
                continue;

            strncpy(charPlaneName, planeName.Data(), 4);
            posVert = TkrClusZ[i];
            posHorAbsExt = ExtrapolateCoordinate(track, posVert);
            invSlope = ( track.GetX2() - track.GetX1() )
                / ( track.GetY2() - track.GetY1() );
            posHorAbs = iView ? TkrClusY[i] : TkrClusX[i];
            // getting the ordinate position from the track in the other view
            const int theOtherView = iView ? 0 : 1;
            if ( refTrack[theOtherView].GetLineStyle() == 1 )
                posHorOrdExt = ExtrapolateCoordinate(refTrack[theOtherView],
                                                     posVert);
            else
                posHorOrdExt = -30;//maybe this should get another default value
            const double res = posHorAbsExt - posHorAbs;
            hRes[iView][iLayer].Fill(res);
            hResiduals[iView].Fill(res);
            resTree.Fill();
        } // end of loop over all clusters
    } // end of loop over all events

    f.Write();
    f.Close();
}

void Residual::DrawResidual(TCut cut) {
    TString canvasName = "DrawResidual";
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c ) {
        //            std::cout << "creating new canvas" << std::endl;
        c = new TCanvas(canvasName, canvasName, 600, 800);
        c->Divide(1,2);
    }

    TFile f(myResFileName);
    TTree* t = (TTree*)f.Get("residualTree");
    t->SetMarkerStyle(20);
    t->SetMarkerSize(0.3);

    c->cd(1);
    gPad->SetTicks(1,1);
    t->Draw("h_abs_ext-h_abs", cut);

    c->cd(2);
    gPad->SetTicks(1,1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    t->Draw("h_abs_ext-h_abs", cut);
    TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->Fit("gaus", "", "", -0.2, 0.2);

    c->cd();
    f.Close();
}

void Residual::DrawResSlope(TCut cut) {
    TString canvasName = "DrawResSlope";
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c ) {
        //            std::cout << "creating new canvas" << std::endl;
        c = new TCanvas(canvasName, canvasName, 600, 800);
        c->Divide(1,2);
    }

    TFile f(myResFileName);
    TTree* t = (TTree*)f.Get("residualTree");
    t->SetMarkerStyle(20);
    t->SetMarkerSize(0.3);

    c->cd(1);
    gPad->SetTicks(1,1);
    t->Draw("h_abs_ext-h_abs:invSlope", cut);

    c->cd(2);
    gPad->SetTicks(1,1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    t->Draw("h_abs_ext-h_abs:invSlope", cut, "prof");
    TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->Fit("pol1", "", "", -0.6, 0.6);

    c->cd();
    f.Close();
}

void Residual::DrawResSlopeAll(TCut cut) {
    TString canvasName = "DrawResSlopeAll";
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c ) {
        c = new TCanvas(canvasName, canvasName, 1100, 800);
        c->Divide(6,6);
    }

    TFile f(myResFileName);
    TTree* t = (TTree*)f.Get("residualTree");
    TMap* myGeometry = myTracker->GetGeometry();
    TMapIter ti(myGeometry);
    TObjString* key;

    gStyle->SetTitleFontSize(0.1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(11);
    //    gStyle->SetStatFontSize(0.1);
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatH(0.25);
    gStyle->SetStatW(0.35);
    gStyle->SetFitFormat(".3f");
    int i=0;
    std::ofstream fout("newGeometry.txt");
    fout.setf(std::ios_base::fixed);
    fout.precision(3);
    while ( (key=(TObjString*)ti.Next()) ) {
        Layer* plane = (Layer*)myGeometry->GetValue(key);
        TString planeName = key->String();
        c->cd(++i);
        gPad->SetTicks(1,1);
        TString onPlane = "name==\"" + planeName + "\"";
        TCut cutOnPlane = onPlane.Data();
        TCut myCut = cut && cutOnPlane;
        t->Draw("h_abs_ext-h_abs:invSlope", myCut, "prof");
        TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
        float dx, dy, dz;
        dx = dy = dz = 0.;
        if ( htemp ) { // histogram might be empty
            htemp->SetTitle(planeName);
            htemp->Fit("pol1", "q", "", -0.6, 0.6);
            gPad->Update();
            TF1* f = htemp->GetFunction("pol1");
            dz = -f->GetParameter(1);
            float dh = f->GetParameter(0);
            if ( plane->GetView() )
                dy = dh;
            else
                dx = dh;
            std::cout << planeName << " shifts: horizonal = " << std::fixed
                      << std::setprecision(3) << dh << " vertical = " << dz
                      << std::endl;
        }
        fout << planeName << ' ' << plane->GetZ() + dz << ' '
             << plane->GetY() + dy << ' ' << plane->GetX() + dx << std::endl;
    }

    fout.close();
    c->cd();
    f.Close();
}

void Residual::DrawResOrd(TCut cut) {
    const TCut cutOrd("h_ord>=0");
    cut += cutOrd;

    TString canvasName = "DrawResOrd";
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c ) {
        c = new TCanvas(canvasName, canvasName, 600, 800);
        c->Divide(1,2);
    }

    TFile f(myResFileName);
    TTree* t = (TTree*)f.Get("residualTree");
    t->SetMarkerStyle(20);
    t->SetMarkerSize(0.3);

    c->cd(1);
    gPad->SetTicks(1,1);
    t->Draw("h_abs_ext-h_abs:h_ord", cut);

    c->cd(2);
    gPad->SetTicks(1,1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    t->Draw("h_abs_ext-h_abs:h_ord", cut, "prof");
    TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->Fit("pol1");

    c->cd();
    f.Close();
}

void Residual::DrawResOrdAll(TCut cut) {
    const TCut cutOrd("h_ord>=0");
    cut += cutOrd;

    TString canvasName = "DrawResOrdAll";
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c ) {
        c = new TCanvas(canvasName, canvasName, 1100, 800);
        c->Divide(6,6);
    }

    TFile f(myResFileName);
    TTree* t = (TTree*)f.Get("residualTree");
    TMap* myGeometry = myTracker->GetGeometry();
    TMapIter ti(myGeometry);
    TObjString* key;

    gStyle->SetTitleFontSize(0.1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(11);
    //    gStyle->SetStatFontSize(0.1);
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatH(0.25);
    gStyle->SetStatW(0.35);
    gStyle->SetFitFormat(".3g");
    int i=0;
    //    std::ofstream fout("newGeometry.txt");
    //    fout.setf(std::ios_base::fixed);
    //    fout.precision(3);
    while ( (key=(TObjString*)ti.Next()) ) {
        Layer* plane = (Layer*)myGeometry->GetValue(key);
        TString planeName = key->String();
        c->cd(++i);
        gPad->SetTicks(1,1);
        TString onPlane = "name==\"" + planeName + "\"";
        TCut cutOnPlane = onPlane.Data();
        TCut myCut = cut && cutOnPlane;
        t->Draw("h_abs_ext-h_abs:h_ord", myCut, "prof");
        TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
        float dx, dy, dz, dangZ;
        dx = dy = dz = 0.;
        if ( htemp ) { // histogram might be empty
            htemp->SetTitle(planeName);
            htemp->Fit("pol1", "q");
            gPad->Update();
            TF1* f = htemp->GetFunction("pol1");
            // dangZ is the rotation around z (check the sign)
            dangZ = f->GetParameter(1);
            // dh is the horizontal shift.  This should have been alredy fitted
            // before.  Most likely, it's 0 or close to 0.
            float dh = f->GetParameter(0);
            if ( plane->GetView() )
                dy = dh;
            else
                dx = dh;
            //   std::cout << planeName << " shifts: horizonal = " << std::fixed
            //             << std::setprecision(3) << dh << " vertical = " << dz
            //                      << std::endl;
        }
        //        fout << planeName << ' ' << plane->GetZ() + dz << ' '
        //      << plane->GetY() + dy << ' ' << plane->GetX() + dx << std::endl;
    }

    //    fout.close();
    c->cd();
    f.Close();
}

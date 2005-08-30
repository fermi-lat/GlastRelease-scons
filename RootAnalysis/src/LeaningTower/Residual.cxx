#include "Residual.h"

ClassImp(Residual)

Residual::Residual(TString filename, TString resFileName, TString geoFileName, int temid)
    : m_temid(temid) {
    if ( geoFileName.Length() == 0 )
        geoFileName = "$ROOTANALYSISROOT/src/LeaningTower/geometry/TowerBgeometry306000517.txt";
    myTracker = new Tracker;
    myTracker->loadGeometry(geoFileName);
    myTracker->SetTower(true);
    myEvent = new Event(filename, myTracker->GetGeometry());
    myResFileName = resFileName;
}

void Residual::Go(int numEntries, int firstEntry) {
    if ( firstEntry < 0 )
        firstEntry = 0;
    const int allEntries = myEvent->GetEntries();
    if ( numEntries == -1 )
        numEntries = allEntries;
    else
        numEntries = std::min(allEntries, firstEntry+numEntries);
    numEntries -= firstEntry;

    std::cout << "Residual::Go(): firstEntry = " << firstEntry
              << " numEntries = " << numEntries << std::endl;

    int bins = 80000;
    float xmin = -400;
    float xmax = 400;
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

    TH1F hResiduals[2];
    hResiduals[0] = TH1F("Xresiduals", "Xresiduals", bins, xmin, xmax);
    hResiduals[1] = TH1F("Yresiduals", "Yresiduals", bins, xmin, xmax);

    TH1F hRes[2][18];
    for ( int i=0; i<18; ++i ) {
        TString title;
        title = TString("X") + (TString("residuals")+=i);
        hRes[0][i] = TH1F(title, title, bins, xmin, xmax);
        title = TString("Y") + (TString("residuals")+=i);
        hRes[1][i] = TH1F(title, title, bins, xmin, xmax);
    }

    const std::vector<TString> planeCols[2] = {
        myTracker->GetPlaneNameCol(0), myTracker->GetPlaneNameCol(1) };

    Recon* recon = myEvent->GetRecon();
    Progress progress;

    // some counters for the statistics
    int usedTrack[2] = { 0, 0 };
    int usedEvent = 0;

    for ( int entry=firstEntry; entry<firstEntry+numEntries; ++entry ) {
        myEvent->Go(entry);

        progress.Go(entry, numEntries, firstEntry);

        if ( myEvent->GetTemId() != m_temid ) 
            continue;

        recon->GetEvent(entry);
        const Int_t TkrNumClus = recon->GetTkrNumClus();
        const Int_t TkrNumTracks = recon->GetTkrNumTracks();
        //        const Int_t TkrTrk1NumClus = recon->GetTkrTrk1NumClus();

        if ( TkrNumTracks != 1 )
            continue;

        // "corrects" cluster positions with respect to the recon results
        int alignSvcFlag = recon->TkrAlignmentSvc(myTracker, true);
        if ( alignSvcFlag )
            continue;

        const Int_t* TkrClusLayer = recon->GetTkrClusLayer();
        const Int_t* TkrClusView = recon->GetTkrClusView();
        const Float_t* TkrClusX = recon->GetTkrClusX();
        const Float_t* TkrClusY = recon->GetTkrClusY();
        const Float_t* TkrClusZ = recon->GetTkrClusZ();

        // helper for the statistics
        bool trackUsed[2] = { false, false };
        bool eventUsed = false;

        // now, loop over all clusters
        for ( int i=0; i<TkrNumClus; ++i ) {
            const int iView = TkrClusView[i];
            const int iLayer = TkrClusLayer[i];
            const TString planeName = GetPlaneName(iLayer, iView);
            std::vector<TString> planeCol = planeCols[iView];
            // remove all hits on tray to which the plane to be studied belongs
            planeCol.erase(find(planeCol.begin(), planeCol.end(), planeName));
            planeCol.erase(find(planeCol.begin(), planeCol.end(),
                                GetTrayTwinPlaneName(iLayer, iView)));

            const TGraph TGclus = recon->GetClusterGraph(planeCol);
            const TLine track = Reconstruct(&TGclus);

            if ( track.GetLineStyle() != 1 ) // no "good" Chi2
                continue;

            strncpy(charPlaneName, planeName.Data(), 4);
            posVert = TkrClusZ[i];
            posHorAbsExt = ExtrapolateCoordinate(track, posVert);
            invSlope = ( track.GetX2() - track.GetX1() )
                / ( track.GetY2() - track.GetY1() );
            posHorAbs = iView ? TkrClusY[i] : TkrClusX[i];
            // getting the ordinate position from the extrapolated position in
            // the other view (filled by Recon::TkrALignmentSvc)
            posHorOrdExt = iView ? TkrClusX[i] : TkrClusY[i];
            const float res = posHorAbsExt - posHorAbs;
            hRes[iView][iLayer].Fill(res);
            hResiduals[iView].Fill(res);
            resTree.Fill();
            trackUsed[iView] = eventUsed = true;
        } // end of loop over all clusters
        if ( trackUsed[0] )
            ++usedTrack[0];
        if ( trackUsed[1] )
            ++usedTrack[1];
        if ( eventUsed )
            ++usedEvent;
    } // end of loop over all events

    std::cout<<"used "<<usedEvent<<" of "<<numEntries<<" events"<<std::endl;
    std::cout << "found " << ++usedTrack[0] << " good tracks in X and "
              << ++usedTrack[1] << " in Y" << std::endl;

    f.Write();
    f.Close();
}

void Residual::DrawResidual(TString plane, TCut cut) {
    const TCut cutPlane("name==\"" + plane + "\"");
    cut += cutPlane;

    TString canvasName = "DrawResidual";
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c ) {
        //            std::cout << "creating new canvas" << std::endl;
        c = new TCanvas(canvasName, canvasName, 600, 800);
        c->Divide(1,2);
    }
 
    gStyle->Reset();
    gStyle->SetOptStat("oeu");
    gStyle->SetOptFit(1111);

    TFile f(myResFileName);
    TTree* t = (TTree*)f.Get("residualTree");
    t->SetMarkerStyle(20);
    t->SetMarkerSize(0.3);

    c->cd(1);
    gPad->SetTicks(1,1);
    gPad->SetLogy();
    t->Draw("h_abs_ext-h_abs");
    TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("residual/mm");

    c->cd(2);
    gPad->SetTicks(1,1);
    t->Draw("h_abs_ext-h_abs", cut);
    htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("residual/mm");
    htemp->Fit("gaus", "", "", -0.2, 0.2);
    gPad->Update();
    TPaveStats* pStats =
        (TPaveStats*)htemp->GetListOfFunctions()->FindObject("stats");
    pStats->SetOptStat(0);
    pStats->SetX1NDC(0.71);
    pStats->SetY1NDC(0.8);
    gPad->Modified();

    c->cd();
    f.Close();
}

void Residual::DrawResSlope(TString plane, TCut cut) {
    const TCut cutPlane("name==\"" + plane + "\"");
    cut += cutPlane;

    TString canvasName = "DrawResSlope";
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c ) {
        //            std::cout << "creating new canvas" << std::endl;
        c = new TCanvas(canvasName, canvasName, 600, 800);
        c->Divide(1,2);
    }

    gStyle->Reset();
    gStyle->SetOptStat("oeu");
    gStyle->SetOptFit(1111);

    TFile f(myResFileName);
    TTree* t = (TTree*)f.Get("residualTree");
    t->SetMarkerStyle(20);
    t->SetMarkerSize(0.3);

    c->cd(1);
    gPad->SetTicks(1,1);
    t->Draw("h_abs_ext-h_abs:invSlope", cut);
    TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("cot(theta)");
    htemp->GetYaxis()->SetTitle("residual/mm");

    c->cd(2);
    gPad->SetTicks(1,1);
    t->Draw("h_abs_ext-h_abs:invSlope", cut, "prof");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("cot(theta)");
    htemp->GetYaxis()->SetTitle("residual/mm");
    htemp->Fit("pol1", "", "", -0.6, 0.6);
    gPad->Update();
    TPaveStats* pStats =
        (TPaveStats*)htemp->GetListOfFunctions()->FindObject("stats");
    pStats->SetOptStat(0);
    pStats->SetX1NDC(0.71);
    pStats->SetY1NDC(0.8);
    gPad->Modified();

    c->cd();
    f.Close();
}

void Residual::DrawResOrd(TString plane, TCut cut) {
    const TCut cutOrd("h_ord>=0");
    const TCut cutPlane("name==\"" + plane + "\"");
    cut += cutOrd + cutPlane;

    TString canvasName = "DrawResOrd";
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c ) {
        c = new TCanvas(canvasName, canvasName, 600, 800);
        c->Divide(1,2);
    }

    gStyle->Reset();
    gStyle->SetOptStat("oeu");
    gStyle->SetOptFit(1111);

    TFile f(myResFileName);
    TTree* t = (TTree*)f.Get("residualTree");
    t->SetMarkerStyle(20);
    t->SetMarkerSize(0.3);

    c->cd(1);
    gPad->SetTicks(1,1);
    t->Draw("h_abs_ext-h_abs:h_ord", cut);
    TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("position in other view/mm");
    htemp->GetYaxis()->SetTitle("residual/mm");

    c->cd(2);
    gPad->SetTicks(1,1);
    t->Draw("h_abs_ext-h_abs:h_ord", cut, "prof");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("position in other view/mm");
    htemp->GetYaxis()->SetTitle("residual/mm");
    htemp->Fit("pol1");
    gPad->Update();
    TPaveStats* pStats =
        (TPaveStats*)htemp->GetListOfFunctions()->FindObject("stats");
    pStats->SetOptStat(0);
    pStats->SetX1NDC(0.71);
    pStats->SetY1NDC(0.8);
    gPad->Modified();

    c->cd();
    f.Close();
}

void Residual::DrawResPos(TString plane, TCut cut) {
    const TCut cutPlane("name==\"" + plane + "\"");
    cut += cutPlane;

    TString canvasName = "DrawResPos";
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c ) {
        c = new TCanvas(canvasName, canvasName, 600, 800);
        c->Divide(1,2);
    }

    gStyle->Reset();
    gStyle->SetOptStat("oeu");
    gStyle->SetOptFit(1111);

    TFile f(myResFileName);
    TTree* t = (TTree*)f.Get("residualTree");
    t->SetMarkerStyle(20);
    t->SetMarkerSize(0.3);

    c->cd(1);
    gPad->SetTicks(1,1);
    t->Draw("h_abs_ext-h_abs:h_abs", cut);
    TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("position/mm");
    htemp->GetYaxis()->SetTitle("residual/mm");

    c->cd(2);
    gPad->SetTicks(1,1);
    t->Draw("h_abs_ext-h_abs:h_abs", cut, "prof");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTitle("position/mm");
    htemp->GetYaxis()->SetTitle("residual/mm");
    htemp->Fit("pol1");
    gPad->Update();
    TPaveStats* pStats =
        (TPaveStats*)htemp->GetListOfFunctions()->FindObject("stats");
    pStats->SetOptStat(0);
    pStats->SetX1NDC(0.71);
    pStats->SetY1NDC(0.8);
    gPad->Modified();

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

    gStyle->Reset();
    gStyle->SetOptFit(11);
    gStyle->SetTitleFontSize(0.1);
    gStyle->SetFitFormat(".3f");
    gStyle->SetStatBorderSize(1);

    TFile f(myResFileName);
    TTree* t = (TTree*)f.Get("residualTree");

    int i=0;
    TList* myGeometry = myTracker->GetGeometry();
    TIter next(myGeometry);
    while ( Layer* plane = (Layer*)next() ) {
        TString planeName = plane->GetName();
        c->cd(++i);
        gPad->SetTicks(1,1);
        TString onPlane = "name==\"" + planeName + "\"";
        TCut cutOnPlane = onPlane.Data();
        TCut myCut = cut && cutOnPlane;
        t->Draw("h_abs_ext-h_abs:invSlope", myCut, "prof");
        TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
        if ( htemp ) { // histogram might be empty
            htemp->GetXaxis()->SetTitle("cot(theta)");
            htemp->GetYaxis()->SetTitle("residual/mm");
            htemp->SetTitle(planeName);
            htemp->Fit("pol1", "q", "", -1.0, 1.0);
            gPad->Update();
            TPaveStats* pStats =
                (TPaveStats*)htemp->GetListOfFunctions()->FindObject("stats");
            pStats->SetOptStat(0);
            pStats->SetX1NDC(0.35);
            pStats->SetY1NDC(0.8);
            gPad->Modified();
            TF1* f = htemp->GetFunction("pol1");
            // if I define the residual as h_abs_ext-h_abs, the horizontal shift
            // is p0, but the vertical is -p1
            const float dz = -f->GetParameter(1);
            plane->SetdZ(damp*dz);
            const float dh = f->GetParameter(0);
            if ( plane->GetView() )
                plane->SetdY(damp*dh);
            else
                plane->SetdX(damp*dh);
            std::cout << planeName << " shifts: horizonal = " << std::fixed
                      << std::setprecision(3) << dh << " vertical = " << dz
                      << std::endl;
        }
    }

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

    gStyle->Reset();
    gStyle->SetOptFit(11);
    gStyle->SetTitleFontSize(0.1);
    gStyle->SetFitFormat(".3g");
    gStyle->SetStatBorderSize(1);

    TFile f(myResFileName);
    TTree* t = (TTree*)f.Get("residualTree");

    int i=0;
    TList* myGeometry = myTracker->GetGeometry();
    TIter next(myGeometry);
    while ( Layer* plane = (Layer*)next() ) {
        TString planeName = plane->GetName();
        c->cd(++i);
        gPad->SetTicks(1,1);
        TString onPlane = "name==\"" + planeName + "\"";
        TCut cutOnPlane = onPlane.Data();
        TCut myCut = cut && cutOnPlane;
        t->Draw("h_abs_ext-h_abs:h_ord", myCut, "prof");
        TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
        if ( htemp ) { // histogram might be empty
            htemp->GetXaxis()->SetTitle("position in other view/mm");
            htemp->GetYaxis()->SetTitle("residual/mm");
            htemp->SetTitle(planeName);
            htemp->Fit("pol1", "q", "", 0, 350);
            gPad->Update();
            TPaveStats* pStats =
                (TPaveStats*)htemp->GetListOfFunctions()->FindObject("stats");
            pStats->SetOptStat(0);
            pStats->SetX1NDC(0.35);
            pStats->SetY1NDC(0.8);
            gPad->Modified();
            TF1* f = htemp->GetFunction("pol1");
            // dangZ is the rotation around z (check the sign)
            const float dangZ = 1000. * f->GetParameter(1);  // mrad
            // dh is the horizontal shift.  This should have been already fitted
            // before.  Most likely, it's 0 or close to 0.
            // However, here it is the offset a hor=0, not the offset for the
            // center of the plane.  Thus, I don't use it here.
            float dh = f->GetParameter(0);
            if ( plane->GetView() )
                plane->SetdRotZ(damp*dangZ);
            else
                plane->SetdRotZ(-damp*dangZ);
            std::cout << planeName << " dh(h=0) = " << std::setprecision(3)
                      << std::fixed << dh << " dRotZ = " << dangZ << std::endl;
        }
    }

    c->cd();
    f.Close();
}

void Residual::DrawResPosAll(TCut cut) {
    std::cout << "MWK: Doesn't do what it claims to!" << std::endl;
    return;
    TString canvasName = "DrawResPosAll";
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c ) {
        c = new TCanvas(canvasName, canvasName, 1100, 800);
        c->Divide(6,6);
    }

    gStyle->Reset();
    gStyle->SetOptFit(11);
    gStyle->SetTitleFontSize(0.1);
    gStyle->SetFitFormat(".3g");
    gStyle->SetStatBorderSize(1);

    TFile f(myResFileName);
    TTree* t = (TTree*)f.Get("residualTree");

    int i=0;
    TList* myGeometry = myTracker->GetGeometry();
    TIter next(myGeometry);
    while ( Layer* plane = (Layer*)next() ) {
        TString planeName = plane->GetName();
        c->cd(++i);
        gPad->SetTicks(1,1);
        TString onPlane = "name==\"" + planeName + "\"";
        TCut cutOnPlane = onPlane.Data();
        TCut myCut = cut && cutOnPlane;
        t->Draw("h_abs_ext-h_abs:h_abs", myCut, "prof");
        TH1F* htemp = (TH1F*)gPad->GetPrimitive("htemp");
        if ( htemp ) { // histogram might be empty
            htemp->GetXaxis()->SetTitle("position/mm");
            htemp->GetYaxis()->SetTitle("residual/mm");
            htemp->SetTitle(planeName);
            htemp->Fit("pol1", "q", "", 0, 350);
            gPad->Update();
            TPaveStats* pStats =
                (TPaveStats*)htemp->GetListOfFunctions()->FindObject("stats");
            pStats->SetOptStat(0);
            pStats->SetX1NDC(0.35);
            pStats->SetY1NDC(0.8);
            gPad->Modified();
            TF1* f = htemp->GetFunction("pol1");
            // dangH is the rotation around a strip (check the sign)
            const float dangH = 1000. * f->GetParameter(1);  // mrad
            // dh is the horizontal shift.  This should have been already fitted
            // before.  Most likely, it's 0 or close to 0.
            // However, here it is the offset a hor=0, not the offset for the
            // center of the plane.  Thus, I don't use it here.
            float dh = f->GetParameter(0);
            // for an x-plane, I determine a rotation around y (and vice versa).
            if ( plane->GetView() )
                plane->SetdRotX(damp*dangH);
            else
                plane->SetdRotY(damp*dangH);
            std::cout << planeName << " dh(h=0) = " << std::setprecision(3)
                      << std::fixed << dh << " dRotH = " << dangH << std::endl;
        }
    }

    c->cd();
    f.Close();
}

void Residual::DrawResOrdCorr(TString plane, TCut cut) {
    const TCut cutOrd("h_ord>=0");
    const TCut cutPlane("name==\"" + plane + "\"");
    cut += cutOrd + cutPlane;

    TString canvasName = "DrawResOrdCorr";
    TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
    if ( !c ) {
        c = new TCanvas(canvasName, canvasName, 600, 800);
        c->Divide(2, 3);
    }

    TFile f(myResFileName);
    TTree* t = (TTree*)f.Get("residualTree");
    t->SetMarkerStyle(20);
    t->SetMarkerSize(0.3);

    TH1F* htemp;

    TString residual;
    TCut cutRes;
    TCut myCut;

    residual = TString("h_abs_ext-h_abs");
    cutRes = TCut("abs(" + residual + ")<0.5");
    myCut = cut + cutRes;
    myCut.Print();

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(11);
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatH(0.25);
    gStyle->SetStatW(0.35);

    c->cd(1);
    gPad->SetTicks(1,1);
    t->Draw(residual+":h_ord", myCut);
    gPad->Update();

    c->cd(3);
    gPad->SetTicks(1,1);
    t->Draw(residual+":h_ord", myCut, "prof");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->Fit("pol1", "q");
    gPad->Update();
    TF1* fun = htemp->GetFunction("pol1");
    float p0 = fun->GetParameter(0);
    float p1 = fun->GetParameter(1);

    c->cd(5);
    gPad->SetTicks(1,1);
    t->Draw(residual, myCut);
    htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->Fit("gaus", "", "", -0.15, 0.15);
    gPad->Update();

    residual += "-(";
    residual += p0;
    residual += ")-h_ord*";
    residual += p1;

    cutRes = TCut("abs(" + residual + ")<0.5");
    myCut = cut + cutRes;
    myCut.Print();

    c->cd(2);
    gPad->SetTicks(1,1);
    t->Draw(residual+":h_ord", myCut);
    gPad->Update();

    c->cd(4);
    gPad->SetTicks(1,1);
    t->Draw(residual+":h_ord", myCut, "prof");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->Fit("pol1", "q", "", 10, 350);
    gPad->Update();

    c->cd(6);
    gPad->SetTicks(1,1);
    t->Draw(residual, myCut);
    htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->Fit("gaus", "q", "", -0.15, 0.15);
    gPad->Update();

    c->cd();
    f.Close();
}

void Residual::SaveGeometry(TString geoFileName) {
    if ( geoFileName == TString() )
        geoFileName = myTracker->GetGeoFileName() + ".new";
    TList* myGeometry = myTracker->GetGeometry();
    std::ofstream fout(geoFileName /*, ios_base::out*/);
    TIter next(myGeometry);
    while ( Layer* plane = (Layer*)next() )
        fout << plane->SaveGeometry() << std::endl;
    fout.close();
}

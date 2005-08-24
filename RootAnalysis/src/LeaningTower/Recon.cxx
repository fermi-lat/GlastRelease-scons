#include "Recon.h"

#include "TF1.h"
#include "TVector3.h"

Recon::Recon(TFile* file,int temid) : m_temid(temid) {
    reconTree = (TTree*)file->Get("Recon");
    if ( reconTree ) {
        reconTree->SetBranchAddress("TkrNumClus",     &TkrNumClus);
        reconTree->SetBranchAddress("TkrNumTracks",   &TkrNumTracks);
        reconTree->SetBranchAddress("TkrTrk1NumClus", &TkrTrk1NumClus);
        reconTree->SetBranchAddress("TkrClusX",        TkrClusX);
        reconTree->SetBranchAddress("TkrClusY",        TkrClusY);
        reconTree->SetBranchAddress("TkrClusZ",        TkrClusZ);
        reconTree->SetBranchAddress("TkrClusLayer",    TkrClusLayer);
        reconTree->SetBranchAddress("TkrClusView",     TkrClusView);
        reconTree->SetBranchAddress("TkrTrk1Clusters", TkrTrk1Clusters);
    }
}

void Recon::GetEvent(int event) {
    int bytes = reconTree->GetEntry(event);
    if ( bytes == 0 )
        std::cerr <<"Recon::GetEvent: couldn't read the Recon part!"<<std::endl;
    // the data are nominally tower 0, which is located at about -700, -700.
    // We would like to have to origin at 0,0 (equivalent to tower 10).
    // There is also a non-understood not-cared-for shift in z.
    //temid is the requested tem to be studied. By default it is bay0

    //Here a systematic translation is applied to bring the tower
    //to 0,0 coordinates (id est bay10). This translation depends on 
    //the position on the grid of the tower under study. 
    //Here we define this translation :
    int offsetX = m_temid%4;
    int offsetY = m_temid/4;
    //    std::cout<<m_temid<<": Translation X and Y "<<offsetX<<" "<<offsetY<<std::endl;
    for ( int i=0; i<TkrNumClus; ++i ) {
        TkrClusX[i] = TkrClusX[i] + 741.09 - offsetX*374.5;
        TkrClusY[i] = TkrClusY[i] + 741.09 - offsetY*374.5;
	//        TkrClusX[i] = TkrClusX[i] + 741.09;
	//        TkrClusY[i] = TkrClusY[i] + 741.09;
        TkrClusZ[i] = TkrClusZ[i] +  18;
    }
}

void Recon::TkrAlignmentSvc(const TList *myGeometry) {
    // this is the "real" alignment.  Z is replaced by the Z of the geometry
    // file.  This seems to be awkward, but the alternative would be to rerun
    // recon every time the position changes.
    // X and Y are taken as corrections, i.e. X and Y of the geometry file are
    // added to the cluster positions.
  
    for ( int i=0; i<TkrNumClus; ++i ) {
        Layer* plane = (Layer*)myGeometry->FindObject(
		          GetPlaneNameFromRecon(TkrClusLayer[i], TkrClusView[i]));
        // translation
        TVector3 v(TkrClusX[i]+plane->GetX(), TkrClusY[i]+plane->GetY(),
                   plane->GetZ());
        /* I had to realize that to make a correction of the rotation I would
           need to know the real 2d in-plane position of the cluster.  However,
           we know the position perpendicular to the strips (i.e. the strip)
           only.  The 2nd coordinate is just the center of the tray.
        // rotation in Z is done with respect to the center of the
        // tower (178.25, 178.25)
        const TVector3 origin(178.25, 178.25, 0);
        v -= origin;
        v.RotateZ(plane->GetRotZ()/1000);
        v += origin;
        */
        // refilling the clusters
        TkrClusX[i] = v.X();
        TkrClusY[i] = v.Y();
        TkrClusZ[i] = v.Z();
    }
}

TGraph Recon::GetAllClustersGraph(const TString view, const int notLayer) const{
    if ( view == "X" )
        return GetAllClustersGraph(0, notLayer);
    else if ( view == "Y" )
        return GetAllClustersGraph(1, notLayer);
    else
        std::cerr << "Recon::GetAllClustersGraph: view = " << view << std::endl;

    std::exit(42);
    return TGraph();
}

TGraph Recon::GetAllClustersGraph(const int view, const int notLayer) const {
    TGraph clusters;
    for ( int i=0; i<GetTkrNumClus(); ++i ) {
        if ( TkrClusView[i] == view && TkrClusLayer[i] != notLayer ) {
            double pos;
            switch ( view ) {
            case 0:
                pos = TkrClusX[i];
                break;
            case 1:
                pos = TkrClusY[i];
                break;
            default:
                pos = 0;
                std::cerr << "Recon::GetAllClustersGraph: view = " << view
                          << std::endl;
                std::exit(42);
            }
            clusters.SetPoint(clusters.GetN(), pos, TkrClusZ[i] );
        }
    }
    return clusters;
}

TGraph Recon::GetClusterGraph(const std::vector<TString> planeCol,
                              bool* exactlyOne) const {
    TGraph tg;
    if ( planeCol.size() == 0 )  // shouldn't happen
        return tg;

    *exactlyOne = true;
    for ( std::vector<TString>::const_iterator it=planeCol.begin();
          it<planeCol.end(); ++it ) {
        // loop over the planes
        int count = 0;
        const int thisView = GetView(*it);
        const int thisReconLayer = GetReconLayer(*it);
        for ( int j=0; j<GetTkrNumClus(); ++j ) {
            // a second loop over all clusters, to select those on the plane
            if ( thisReconLayer!=TkrClusLayer[j] || thisView!=TkrClusView[j] )
                continue; // cluster is not on this plane
            ++count;
            double pos = 0;
            if ( thisView == 0 )
                pos = TkrClusX[j];
            else if ( thisView == 1 )
                pos = TkrClusY[j];
            tg.SetPoint(tg.GetN(), pos, TkrClusZ[j] );
        }
        if ( count != 1 )
            *exactlyOne = false;
    }
    return tg;
}

TGraph Recon::GetTrk1ClustersGraph(const TString view, const int notLayer)const{
    if ( view == "X" )
        return GetTrk1ClustersGraph(0, notLayer);
    else if ( view == "Y" )
        return GetTrk1ClustersGraph(1, notLayer);
    else
        std::cerr << "Recon::GetTrk1ClustersGraph: view = " << view <<std::endl;

    std::exit(42);
    return TGraph();
}

TGraph Recon::GetTrk1ClustersGraph(const int view, const int notLayer) const {
    TGraph clusters;
    std::cout<<"test: view and notlayer requested"<< view<<" "<<notLayer<<std::endl;
    std::cout<<"test: GetTkrTrk1NumClus"<<GetTkrTrk1NumClus()<<std::endl;
    for ( int j=0; j<GetTkrTrk1NumClus(); ++j ) {
        int i = TkrTrk1Clusters[j];
	std::cout<<"test: j i view layer"<<j<<" "<<i<<" "<<TkrClusView[i]<<" "<<TkrClusLayer[i]<<std::endl;
        if ( TkrClusView[i] == view && TkrClusLayer[i] != notLayer ) {
            double pos;
            switch ( view ) {
            case 0:
                pos = TkrClusX[i];
                break;
            case 1:
                pos = TkrClusY[i];
                break;
            default:
                pos = 0;
                std::cerr << "Recon::GetTrk1ClustersGraph: view = " << view
                          << std::endl;
                std::exit(42);
            }
            clusters.SetPoint(clusters.GetN(), pos, TkrClusZ[i] );
        }
    }
    return clusters;
}


TLine Reconstruct(const TGraph* XY) {
    // LineStyle 0: less than two points
    // LineStyle 1: fit with good chi2
    // LineStyle 2: fit with bad chi2

    int N = XY->GetN();
    if ( N < 2 ) {
        TLine invalid;
        invalid.SetLineStyle(0);
        return invalid;
    }

    double* X = XY->GetX();
    double* Y = XY->GetY();
    double MinChi2 = 1E6;
    double A = 0;
    double B = 0;

    // comment the code that improves the fit but eats CPU
    //    if ( N <= 3 ) {
        // swapping columns.  This improves the fit!
        TGraph g(N, Y, X);
        g.Fit("pol1", "Q");
        TF1* f = (TF1*)g.FindObject("pol1");
        MinChi2 = f->GetChisquare();      
        A = f->GetParameter(0);
        B = f->GetParameter(1);
        delete f;
        /*
    }
    else {
        for ( int i=0; i<N ; ++i ) {
            // swapping columns.  This improves the fit!
            TGraph g(N, Y, X);
            g.RemovePoint(i);
            g.Fit("pol1", "Q");
            TF1* f = (TF1*)g.FindObject("pol1");
            const double chi2 = f->GetChisquare();
            if ( chi2 < MinChi2 ) {
                MinChi2 = chi2;
                A = f->GetParameter(0);
                B = f->GetParameter(1);
            }
            delete f;
        }
    }
        */

    const double Y1 = 0.0;
    const double Y2 = 700.0;
    // remember: x and y are swapped!
    const double x1 = A + B * Y1;
    const double y1 = Y1;
    const double x2 = A + B * Y2;
    const double y2 = Y2;
    //    std::cout << x1 << ' ' << y1 << ' ' << x2 << ' ' << y2 << std::endl;

    TLine track(x1, y1, x2, y2);
    track.SetLineColor(2);

    static const double MaxChi2 = 1;
    //    std::cout << "chisqr " << MinChi2 << std::endl;
    if ( MinChi2 > MaxChi2 )
        track.SetLineStyle(2);

    return track;
}

bool IsValid(const TLine& l) { return l.GetLineStyle(); }

double ExtrapolateCoordinate(const TLine& track, const double z) {
    const double x1 = track.GetX1();
    const double x2 = track.GetX2();
    const double y1 = track.GetY1();
    const double y2 = track.GetY2();
    const double newX = ( x2 == x1 )  ?  x1  :  x1 + (z-y1) / (y2-y1) * (x2-x1);
    return newX;
}

ClassImp(Recon)

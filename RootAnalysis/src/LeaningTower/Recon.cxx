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
    for ( int i=0; i<TkrNumClus; ++i ) {
        TkrClusX[i] = TkrClusX[i] + ( 2 - m_temid % 4 ) * TOWERPITCH;
        TkrClusY[i] = TkrClusY[i] + ( 2 - m_temid / 4 ) * TOWERPITCH;
        TkrClusZ[i] = TkrClusZ[i] +  18;
    }
}

int Recon::TkrAlignmentSvc(const Tracker *myTracker, const bool rotate) {
    // this is the "real" alignment.  Z is replaced by the Z of the geometry
    // file.  This seems to be awkward, but the alternative would be to rerun
    // recon every time the position changes.
    // X and Y are taken as corrections, i.e. X and Y of the geometry file are
    // added to the cluster positions.

    int alignSvcFlag = 0;

    const TList* myGeometry = myTracker->GetGeometry();

    // translation
    for ( int i=0; i<TkrNumClus; ++i ) {
        Layer* plane = (Layer*)myGeometry->FindObject(
                        GetPlaneNameFromRecon(TkrClusLayer[i], TkrClusView[i]));
        /*
        TVector3 v(TkrClusXorig[i]+plane->GetX(), TkrClusYorig[i]+plane->GetY(),
                   plane->GetZ());
        TkrClusX[i] = v.X();
        TkrClusY[i] = v.Y();
        TkrClusZ[i] = v.Z();
        */
        TkrClusX[i] = TkrClusXorig[i] + plane->GetX();
        TkrClusY[i] = TkrClusYorig[i] + plane->GetY();
        TkrClusZ[i] = plane->GetZ();
    }

    // rotation
    if ( rotate ) {
        // let's first define two tracks for later determination of the
        // approximate position in the plane of interest.
        // We don't care that there is exactly one cluster per plane.  If it
        // would matter, the Chi2 would be bad anyway!
        TLine refTrack[2];
        for ( int i=0; i<2; ++i ) {
            // adding all clusters on planeCols[i] planes to TGclus
            const TGraph TGclus =GetClusterGraph(myTracker->GetPlaneNameCol(i));
            refTrack[i] = Reconstruct(&TGclus);
            // If the line style is != 1, the ChiSqr of the track is worse than
            // 1.  Let's flag this in the return value of the alignment service.
            if ( refTrack[i].GetLineStyle() != 1 )
                alignSvcFlag += ( 1 << i );
        }

        for ( int i=0; i<TkrNumClus; ++i ) {
            Layer* plane = (Layer*)myGeometry->FindObject(
                               GetPlaneNameFromRecon(TkrClusLayer[i],
                                                     TkrClusView[i]));
            TVector3 v(TkrClusX[i], TkrClusY[i], TkrClusZ[i]);
            const int theOtherView = TkrClusView[i] ? 0 : 1;
            // let's save the extrapolated position as the cluster's real
            // position (as opposed to the planes center).
            v[theOtherView] = ExtrapolateCoordinate(refTrack[theOtherView],
                                                    TkrClusZ[i]);
            float rot;
            rot = plane->GetRotZ();
            //            v.Print();
            if ( rot ) {
                //                std::cout << "rotation z " << rot << std::endl;
                const TVector3 axisOffset(178.25, 178.25, v[2]);
                v -= axisOffset;
                v.RotateZ(rot/1000);
                v += axisOffset;
                //                v.Print();
            }

            /* MWK: Does yet do what it claims to do!
            rot = plane->GetRotY();
            if ( rot ) {
                //                std::cout << "rotation y " << rot << std::endl;
                const TVector3 axisOffset(178.25, v[1], v[2]);
                v -= axisOffset;
                v.RotateY(rot/1000);
                v += axisOffset;
                //                v.Print();
            }

            rot = plane->GetRotX();
            if ( rot ) {
                //                std::cout << "rotation x " << rot << std::endl;
                const TVector3 axisOffset(v[0], 178.25, v[2]);
                v -= axisOffset;
                v.RotateX(rot/1000);
                v += axisOffset;
                //                v.Print();
                //                std::exit(42);
            }
            */

            // refilling the clusters
            TkrClusX[i] = v.X();
            TkrClusY[i] = v.Y();
            TkrClusZ[i] = v.Z();
        }
    }
    return alignSvcFlag;
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
            float pos;
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
            float pos = 0;
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
    //    std::cout<<"test: view and notlayer requested"<< view<<" "<<notLayer<<std::endl;
    //    std::cout<<"test: GetTkrTrk1NumClus"<<GetTkrTrk1NumClus()<<std::endl;
    for ( int j=0; j<GetTkrTrk1NumClus(); ++j ) {
        int i = TkrTrk1Clusters[j];
        //	std::cout<<"test: j i view layer"<<j<<" "<<i<<" "<<TkrClusView[i]<<" "<<TkrClusLayer[i]<<std::endl;
        if ( TkrClusView[i] == view && TkrClusLayer[i] != notLayer ) {
            float pos;
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

    Double_t* X = XY->GetX();
    Double_t* Y = XY->GetY();
    Double_t MinChi2 = 1E6;
    Double_t A = 0;
    Double_t B = 0;

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
            const Double_t chi2 = f->GetChisquare();
            if ( chi2 < MinChi2 ) {
                MinChi2 = chi2;
                A = f->GetParameter(0);
                B = f->GetParameter(1);
            }
            delete f;
        }
    }
        */

    const Double_t Y1 = 0.0;
    const Double_t Y2 = 700.0;
    // remember: x and y are swapped!
    const Double_t x1 = A + B * Y1;
    const Double_t y1 = Y1;
    const Double_t x2 = A + B * Y2;
    const Double_t y2 = Y2;
    //    std::cout << x1 << ' ' << y1 << ' ' << x2 << ' ' << y2 << std::endl;

    TLine track(x1, y1, x2, y2);
    track.SetLineColor(2);

    static const Double_t MaxChi2 = 1;
    //    std::cout << "chisqr " << MinChi2 << std::endl;
    if ( MinChi2 > MaxChi2 )
        track.SetLineStyle(2);

    return track;
}

bool IsValid(const TLine& l) { return l.GetLineStyle(); }

float ExtrapolateCoordinate(const TLine& track, const float z) {
    const Double_t x1 = track.GetX1();
    const Double_t x2 = track.GetX2();
    const Double_t y1 = track.GetY1();
    const Double_t y2 = track.GetY2();
    const Double_t newX = ( x2 == x1 ) ? x1 : x1 + (z-y1) / (y2-y1) * (x2-x1);
    return (float)newX;
}

ClassImp(Recon)

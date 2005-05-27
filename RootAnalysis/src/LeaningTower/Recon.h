#ifndef __LEANINGTOWER_RECON__
#define __LEANINGTOWER_RECON__

#include "Layer.h"

#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLine.h"
#include "TList.h"

#include <iomanip>
#include <iostream>
#include <vector>

class Recon : public TObject {
 public:
    Recon(TFile*);

    void GetEvent(int);
    void TkrAlignmentSvc(const TList* myGeometry);
    Int_t GetEntries()        const { return (Int_t)reconTree->GetEntries(); }
    Int_t GetTkrNumClus()     const { return TkrNumClus; }
    Int_t GetTkrNumTracks()   const { return TkrNumTracks; }
    Int_t GetTkrTrk1NumClus() const { return TkrTrk1NumClus; }
    const Int_t* GetTkrClusLayer() const { return TkrClusLayer; }
    const Int_t* GetTkrClusView()  const { return TkrClusView; }
    const Float_t* GetTkrClusX()  const { return TkrClusX; }
    const Float_t* GetTkrClusY()  const { return TkrClusY; }
    const Float_t* GetTkrClusZ()  const { return TkrClusZ; }
    const Int_t* GetTkrTrk1Clus()  const { return TkrTrk1Clusters; }
    TGraph GetAllClustersGraph(const TString view, const int notLayer=-1) const;
    TGraph GetAllClustersGraph(const int view, const int notLayer=-1) const;
    TGraph GetClusterGraph(const std::vector<TString> planeCol,
                           bool* exactlyOne) const;
    TGraph GetClusterGraph(const std::vector<TString> planeCol) const {
        bool dummy;
        return GetClusterGraph(planeCol, &dummy);
    }
    TGraph GetClusterGraph(const TString plane) const {
        std::vector<TString> planeCol;
        planeCol.push_back(plane);
        return GetClusterGraph(planeCol);
    }
    TGraph GetTrk1ClustersGraph(const TString view, const int notLayer=-1)const;
    TGraph GetTrk1ClustersGraph(const int view, const int notLayer=-1) const;

    void PrintTkrCluster(int i) const {
        std::cout << std::setw(3) << i << std::setw(3) << TkrClusLayer[i]
                  << std::setw(2) << TkrClusView[i] << ' ' << TkrClusX[i] << ' '
                  << TkrClusY[i] << ' ' << TkrClusZ[i] << std::endl;
    }
    void PrintTkrClusters() const {
        for ( int i=0; i<GetTkrNumClus(); ++i )
            PrintTkrCluster(i);
    }
    void PrintTkrTrk1Clusters() const {
        for ( int i=0; i<GetTkrTrk1NumClus(); ++i )
            PrintTkrCluster(TkrTrk1Clusters[i]);
    }

 private:
    TTree* reconTree;

    // attention: the maximum length of the arrays can be 128, but in the tree
    // are being stored exactly TkrNumClus/TkrTrk1NumClus.

    // clusters
    Int_t TkrNumClus;
    Int_t TkrClusLayer[128];
    Int_t TkrClusView[128];
    Float_t TkrClusX[128];
    Float_t TkrClusY[128];
    Float_t TkrClusZ[128];
    // tracks
    Int_t TkrNumTracks;
    Int_t TkrTrk1NumClus;
    Int_t TkrTrk1Clusters[128];

    // last line of class def
    ClassDef(Recon, 1)
};

TLine Reconstruct(const TGraph*);
bool IsValid(const TLine& l);
//bool IsInvalid(const TLine& l) { return l.GetX1() == 0.0 && l.GetX2() == 0.0
//                                    && l.GetY1() == 0.0 && l.GetY2() == 0.0; }
double ExtrapolateCoordinate(const TLine&, const double z);

#endif

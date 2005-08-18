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
    bool debug;
    int m_temid;
public:
    Efficiency(const TString="MyRootFile.root", const TString="efficiency.root",
               TString geoFileName="", bool debug=false, int temid=0);
    virtual ~Efficiency() {
        delete myEvent;
        delete myTracker;
    }
    // structure to save the ntuple values before filling the tree
    struct Ntuple {
        Int_t eventId;
        Char_t charPlaneName[4];
        Float_t xExt, yExt, siDist, res;
        Float_t siDistOther, resOther;
        Bool_t trigger;
    };

    void Go(int lastEntry=-1);
    // if Get_3_in_a_row would work on the array of TriggerDiagnostics in the
    // Header tree, it would be faster.  But, currently, this array is in the
    // order of GTCC/GTRC.  This stuff should get extracted from TreeMaker
    // and here, and go into it's own class.
    bool Get_3_in_a_row(const int layer) const;
    void setDebug(const bool d) { debug = d; }

    void Draw2D(const TString planeName="all", TCut="",
                const float residualDist=1, const float borderWidth=1) const;
    void GetEfficiency(const TString planeName="all", const TCut="",
                   const float residualDist=1, const float borderWidth=1) const;
    void DrawEfficiency(const TString planeName, TCut="",
                   const float residualDist=1, const float borderWidth=1) const;
private:
    void PrintEfficiency(TString planeName, const int hits,
                         const int missing, const int all=0) const;
    // I use functions to define cuts.  This enforces consistency over the
    // various methods.
    // elements
    TCut isPlane(const TString planeName) const {
        return TCut(TString("plane==\"") + planeName + TString("\""));
    }
    // position
    TCut inActiveArea() const { return TCut("siDist<0"); }
    TCut deepInActiveArea(const float borderWidth) const {
        // not only in active area, but also not close to the edge!
        return TCut(TString("siDist<-")+=borderWidth);
    }
    TCut closeToEdge(const float borderWidth) const {
        return inActiveArea() + !deepInActiveArea(borderWidth);
    }
    // hit found
    TCut hitFoundInThisView(const float residualDist) const {
        return TCut(TString("res<")+=residualDist);
    }
    TCut hitFoundInOtherView(const float residualDist) const {
        return TCut(TString("resOther<")+=residualDist);
    }
    // missing hit definition
    TCut missingHit(const float residualDist) const {
        // missingHit is not equal to !hitFoundInThisView !!!!
        // a hit is missing, if:
        //   1) we don't find a hit in the plane within residualDist
        //   2) we find a hit in the other plane of the same layer within
        //      residualDist.
        // If also in the other plane the hit is missing, it is very probable
        // that the particle was stopped.  But, we continue tracing the track
        // through the tower.  Hopefully, we don't get screwed by noise hits.
        return !hitFoundInThisView(residualDist)
            && hitFoundInOtherView(residualDist);
    }
    // ********** compound cuts **********
    // hit found deep in active area (not close to the edge)
    TCut hitFoundDeepInActiveArea(const float residualDist,
                                  const float borderWidth) const {
        return deepInActiveArea(borderWidth)
            && hitFoundInThisView(residualDist);
    }
    // hit found close to the edge (in active area)
    TCut hitFoundCloseToEdge(const float residualDist, const float borderWidth)
        const {
        return closeToEdge(borderWidth) && hitFoundInThisView(residualDist);
    }
    // missing hit deep in active area (not close to the edge)
    TCut missingHitDeepInActiveArea(const float residualDist,
                                    const float borderWidth) const {
        return deepInActiveArea(borderWidth)
            && missingHit(residualDist);
    }
    // missing hit close to the edge (in active area)
    TCut missingHitCloseToEdge(const float residualDist,
                               const float borderWidth) const {
        return closeToEdge(borderWidth) && missingHit(residualDist);
    }
    // rejected missing hits
    TCut rejectedMissingHit(const float residualDist) const {
        return inActiveArea() && !hitFoundInThisView(residualDist)
            && !missingHit(residualDist);
    }
    
    ClassDef(Efficiency,1)
};


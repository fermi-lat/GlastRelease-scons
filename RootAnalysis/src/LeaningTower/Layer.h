#ifndef __LEANINGTOWER_LAYER__
#define __LEANINGTOWER_LAYER__

#include "TLine.h"
#include "TText.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"

#include <vector>
#include <iostream>

//////// HISTORICAL //////////
// What is called "Layer" here, is in Ritz notation a plane!
// A plane (here "Layer") is a single silicon plane.
// We have 36 of these.
//////// HISTORICAL //////////

class Layer : public TObject {
 public:
    Layer(TString name, double height);
    ~Layer();
    void SetPlanesForFittingCol(std::vector<TString> v) { planesForFittingCol = v; }

    TString GetLayerName() const { return Name; }
    double GetHeight()     const { return Z; }

    double GetCoordinate(int stripId);
    bool checkActiveArea(double ParallelCoordinate, double NormalCoordinate, double BorderWidth);
    void DrawLayer();
    void DrawGhostLayer();
    //  bool IsX() const { return Name.Contains("X"); }
    //  bool IsY() const { return Name.Contains("Y"); }
    // the following 3 definitions rely on names of the structure "X7" or "Y13"
    bool IsX() const { return Name(0,1) == "X"; }
    bool IsY() const { return Name(0,1) == "Y"; }
    int GetLayer() const { return atoi(Name(1,Name.Length()-1).Data()); }
    int GetReconLayer() const { return 17 - GetLayer(); }  // digi counts up, recon down
    int GetView() const {
        if ( IsX() ) return 0;
        if ( IsY() ) return 1;
        return -1;
    }

    TString Name;
    double ShiftX;
    double ShiftY;
    double Z;
    TLine *LayerLine;
    TLine *LadderLine[4];
    TTree *LayerTree;
    TText *LayerLabel;
    void SetTree(TFile *file);
    void GetEvent(int event);

    int TkrNumHits;
    int TkrHits[128];

    void AddHitInActiveArea() { HitsInActiveArea++; }
    void AddMissedHit() { MissedHits++; }
    double GetInefficiency() {
        if(HitsInActiveArea)
            return (double) MissedHits/HitsInActiveArea;
    }
    double GetEfficiency() { return 1.0-GetInefficiency(); }
    Bool_t GetTriggerReq(Bool_t side) { return side ? TriggerReq1 : TriggerReq0; }
    const std::vector<TString>& GetPlanesForFittingCol() const { return planesForFittingCol; }
  
 private:
    double EDGE_WIDTH,WAFER_WIDTH,STRIP_PITCH,LADDER_SEPARATION;
    int HitsInActiveArea,MissedHits;
    Bool_t TriggerReq0;
    Bool_t TriggerReq1;
    std::vector<TString> planesForFittingCol;

    // last line
    ClassDef(Layer, 1)
};

// functions which deal with layers, but who don't need an layer object

inline int GetView(const TString name) {
    if ( name(0,1) == "X" )
        return 0;
    if ( name(0,1) == "Y" )
        return 1;
    return -1;
}
inline int GetLayer(const TString name) { return atoi(name(1,name.Length()-1).Data()); }
inline int GetReconLayer(const TString name) { return 17 - GetLayer(name); }

inline int GetTray(const int layer, const int view) { return layer + ( layer + view + 1 ) % 2; }
inline int GetTray(const TString name) { return GetTray(GetLayer(name), GetView(name)); }

TString GetPlaneName(const int layer, const int view) {
    TString v = view ? "Y" : "X";
    return v += layer;
}

TString GetTwinPlaneName(const int layer, const int view) {
    const int tray = GetTray(layer, view);
    // tray with index i: bottom plane is on layer i-1, top on layer i
    const int newLayer = ( tray == layer ) ? layer - 1 : layer + 1;
    if ( newLayer < 0 || newLayer > 17 )
        return "";
    return GetPlaneName(newLayer, view);
}
TString GetTwinPlaneName(const TString name) { return GetTwinPlaneName(GetLayer(name), GetView(name)); }

#endif

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

const float TOWERPITCH        = 374.5;
const float SIWAFERSIDE       = 89.5;
const float SIWAFERACTIVESIDE = 87.552;
const float STRIPPITCH        =  0.228;
const float LADDERGAP         =  0.2; 
const float SSDGAP            =  0.025;

const float INACTIVEBORDERWIDTH = 0.5 * ( SIWAFERSIDE - SIWAFERACTIVESIDE );
const float TOWERPITCH2         = 0.5 * TOWERPITCH;
const float SIWAFERSIDE2        = 0.5 * SIWAFERSIDE;
const float SIWAFERACTIVESIDE2  = 0.5 * SIWAFERACTIVESIDE;

class Layer : public TNamed {
 public:
    Layer(TString name, float pz, float py=0, float px=0,
          float rz=0, float ry=0, float rx=0);
    ~Layer();
    void SetPlanesForFittingCol(std::vector<TString> v) {planesForFittingCol=v;}

    float GetHeight() const { return GetZ(); }
    float GetZ()      const { return Z; }
    float GetY()      const { return Y; }
    float GetX()      const { return X; }
    float GetRotZ()   const { return rotZ; }
    float GetRotY()   const { return rotY; }
    float GetRotX()   const { return rotX; }
    void SetdZ(float a)     { dZ = a; }
    void SetdY(float a)     { dY = a; }
    void SetdX(float a)     { dX = a; }
    void SetdRotZ(float a)  { dRotZ = a; }
    void SetdRotY(float a)  { dRotY = a; }
    void SetdRotX(float a)  { dRotX = a; }
    std::string SaveGeometry() const;

    float GetTemX(int i)       const { return ( i % 4 - 2 ) * TOWERPITCH; }
    float GetTemY(int i)       const { return ( i / 4 - 2 ) * TOWERPITCH; }
    float GetTemCenX(int i)    const { return GetTemX(i) + TOWERPITCH2; }
    float GetTemCenY(int i)    const { return GetTemY(i) + TOWERPITCH2; }
    // ladder and wafer coordinates are LOCAL to the plane!
    float GetLadderCenX(int i) const {
        return TOWERPITCH2 + ( i - 1.5 ) * ( SIWAFERSIDE + LADDERGAP ); }
    float GetWaferCenY(int i)  const {
        return TOWERPITCH2 + ( i - 1.5 ) * ( SIWAFERSIDE + SSDGAP ); }
    float GetLadderXmin(int i) const { return GetLadderCenX(i) - SIWAFERSIDE2; }
    float GetLadderXmax(int i) const { return GetLadderCenX(i) + SIWAFERSIDE2; }
    float GetWaferYmin(int i)  const { return GetWaferCenY(i) - SIWAFERSIDE2; }
    float GetWaferYmax(int i)  const { return GetWaferCenY(i) + SIWAFERSIDE2; }

    float GetCoordinate(int stripId);
    float activeAreaDist(const float x, const float y) const;
    void DrawLayer();
    void DrawGhostLayer();
    // the following 3 definitions rely on names of the structure "X7" or "Y13"
    bool IsX() const { return fName(0,1) == "X"; }
    bool IsY() const { return fName(0,1) == "Y"; }
    int GetLayer() const { return atoi(fName(1,fName.Length()-1).Data()); }
    // digi counts up, recon down
    int GetReconLayer() const { return GetLayer(); }
    int GetView() const {
        if ( IsX() ) return 0;
        if ( IsY() ) return 1;
        return -1;
    }
    void SetTree(TFile *file);
    void GetEvent(int event) { LayerTree->GetEntry(event); }
    Bool_t GetTriggerReq(Bool_t side) { return side ? TriggerReq1 :TriggerReq0;}
    const std::vector<TString>& GetPlanesForFittingCol() const {
        return planesForFittingCol; }

    friend class Event;

 private:
    float X, Y, Z, dX, dY, dZ;
    float rotX, rotY, rotZ, dRotX, dRotY, dRotZ;
    TLine *LayerLine;
    TLine *LadderLine[4];
    TTree *LayerTree;
    TText *LayerLabel;

    int TkrNumHits;
    int TkrHits[128];
    int ToT0, ToT1;

    Bool_t TriggerReq0;
    Bool_t TriggerReq1;
    std::vector<TString> planesForFittingCol;

    // last line of class def
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
inline int GetLayer(const TString s) { return atoi(s(1,s.Length()-1).Data()); }
//JCT inline int GetReconLayer(const TString name) { return GetLayer(name); }
inline int GetReconLayer(const TString name) { return GetLayer(name); }

inline int GetTray(const int l, const int v) { return l + ( l + v + 1 ) % 2; }
inline int GetTray(const TString s) { return GetTray(GetLayer(s), GetView(s)); }

inline TString GetPlaneName(const int layer, const int view) {
    TString v = view ? "Y" : "X";
    return v += layer;
}

inline TString GetPlaneNameFromRecon(const int l, const int v) {
    return GetPlaneName(l, v); }

/// finds the name of the other plane on the same tray
inline TString GetTrayTwinPlaneName(const int layer, const int view) {
    const int tray = GetTray(layer, view);
    // tray with index i: bottom plane is on layer i-1, top on layer i
    const int newLayer = ( tray == layer ) ? layer - 1 : layer + 1;
    if ( newLayer < 0 || newLayer > 17 )
        return "";
    return GetPlaneName(newLayer, view);
}
inline TString GetTrayTwinPlaneName(const TString name) {
    return GetTrayTwinPlaneName(GetLayer(name), GetView(name)); }

/// finds the name of the other plane on the same layer
inline TString GetLayerTwinPlaneName(const int layer, const int view) {
 return GetPlaneName(layer, view?0:1); }
inline TString GetLayerTwinPlaneName(const TString name) {
    return GetLayerTwinPlaneName(GetLayer(name), GetView(name)); }

#endif

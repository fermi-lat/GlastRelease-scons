#ifndef __LEANINGTOWER_TRACKER__
#define __LEANINGTOWER_TRACKER__

#include "TMap.h"
#include "TString.h"
#include "TCanvas.h"

#include <vector>

class Tracker {
 public:
    Tracker();
    virtual ~Tracker();

    void loadGeometry(TString filename="src/LeaningTower/geometry/stack2geometry.txt");
    void loadFitting(TString filename="src/LeaningTower/geometry/FittingPlanes.txt");
    void Display(TCanvas*);
    TMap* GetGeometry() const { return myGeometry; }
    std::vector<TString> GetPlaneNameCol(const int view, const bool onlyTheGood=false) const;
    std::vector<TString> GetPlaneNameCol(const TString view, const bool onlyTheGood=false) const;
    void IsTower(bool tower) {TOWER=tower;}
 private:
    bool TOWER;
    TMap* myGeometry;

    // last line
    ClassDef(Tracker, 1)
};

#endif

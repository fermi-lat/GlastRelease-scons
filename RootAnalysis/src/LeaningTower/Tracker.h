#ifndef __LEANINGTOWER_TRACKER__
#define __LEANINGTOWER_TRACKER__

#include "TList.h"
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"

#include <vector>

class Tracker {
 public:
    Tracker();
    virtual ~Tracker();
    void loadGeometry(TString=
              "$ROOTANALYSISROOT/src/LeaningTower/geometry/Tower0Geometry.txt");
    void loadFitting(TString=
               "$ROOTANALYSISROOT/src/LeaningTower/geometry/FittingPlanes.txt");
    void Display(TCanvas*);
    TList* GetGeometry() const { return myGeometry; }
    std::vector<TString> GetPlaneNameCol(const int view,
                                         const bool onlyTheGood=false) const;
    std::vector<TString> GetPlaneNameCol(const TString view,
                                         const bool onlyTheGood=false) const;
    void SetTower(const bool tower) { TOWER = tower; }
 private:
    bool TOWER;
    TList* myGeometry;

    // last line of class def
    ClassDef(Tracker, 1)
};

#endif

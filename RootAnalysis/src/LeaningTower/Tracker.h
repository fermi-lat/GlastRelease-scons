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
    void loadGeometry(TString);
    void Display(TCanvas*);
    TList* GetGeometry() const { return myGeometry; }
    TString GetGeoFileName() const { return m_geoFileName; }
    std::vector<TString> GetPlaneNameCol(const int view) const;
    std::vector<TString> GetPlaneNameCol(const TString view) const;
    void SetTower(const bool tower) { TOWER = tower; }
 private:
    bool TOWER;
    TList* myGeometry;
    TString m_geoFileName;

    // last line of class def
    ClassDef(Tracker, 1)
};

#endif

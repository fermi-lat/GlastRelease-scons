#include "TCut.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLine.h"
#include "TList.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "Tracker.h"
#include "Layer.h"
#include "Event.h"
#include "Recon.h"
#include "Progress.h"

const float damp = 0.63;

class Residual {
 private:
    Event*   myEvent;
    Tracker* myTracker;
    TString  myResFileName;
    int m_temid;
 public:
    Residual(const TString="MyRootFile.root", const TString="residual.root",
             const TString geo="", int temid=0);
    virtual ~Residual() {
        delete myEvent;
        delete myTracker;
    }

    void Go(int numEntries=-1, int firstEntry=0);
    // general for DrawXxx: residual is defined as h_abs_ext-h_abs
    // ***********
    // draws the residual, top without fitting, bottom with gauss fit
    void DrawResidual(TString plane, TCut="");
    // draws top the residual vs. inv. slope, bottom the profile with pol1 fit
    void DrawResSlope(TString plane, TCut="");
    // draws top the residual vs. other coordinate, bottom profile with pol1 fit
    void DrawResOrd(TString plane, TCut="");
    // draws top the residual vs. position, bottom profile with pol1 fit
    void DrawResPos(TString plane, TCut="");
    // draws the slope profiles for all planes
    void DrawResSlopeAll(TCut="abs(h_abs_ext-h_abs)<1");
    // draws the ordinate profiles for all planes
    void DrawResOrdAll(TCut="abs(h_abs_ext-h_abs)<1");
    // draws the position profiles for all planes
    void DrawResPosAll(TCut="abs(h_abs_ext-h_abs)<1");
    // attempt to do a "statistical" (not event be event) correction of rotZ of
    // planes.  2x3 histograms.  Left uncorrected, right rotZ corrected.  Top
    // ordinate vs. residual, middle profile of that, bottom residual.
    void DrawResOrdCorr(TString plane, TCut="");
    // saves the geometry into geoFileName
    void Residual::SaveGeometry(TString geoFileName="");

    ClassDef(Residual,1)
};

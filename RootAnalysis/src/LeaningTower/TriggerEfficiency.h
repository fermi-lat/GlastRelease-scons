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
#include "Progress.h"

#include <bitset>

// I simplified the code a lot using the new GetTkrDigi3RowBits() and
// GetTkrTrgReq3RowBits() functions.  Thus, I could remove a lot code.  However,
// to log the deviation of TkrDigis and TriggerReqs, I need to keep some of it.
// I you want to check for differences, define TRIGGER_EFF_DEBUG.
#define TRIGGER_EFF_DEBUG

class TriggerEfficiency {
private:
    Event*   myEvent;
    Tracker* myTracker;
public:
    TriggerEfficiency(const TString="MyRootFile.root", TString geoFileName="");
    virtual ~TriggerEfficiency() {
        delete myEvent;
        delete myTracker;
    }

    void Go(int lastEntry=-1);

    ClassDef(TriggerEfficiency,1)
};

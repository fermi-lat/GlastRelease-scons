#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TList.h"
#include "TObjectTable.h"
#include "TObjString.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"

#include "Tracker.h"
#include "Layer.h"
#include "Event.h"
#include "Recon.h"

#include <algorithm>

class EventDisplay {
private:
    Event*   myEvent;
    Tracker* myTracker;
    TCanvas* myEventDisplay;

    TText LabelNumHits[36];
    TText TriggerReqText[36][2];

    TGraph *TGraphXhits;
    TGraph *TGraphYhits;

    TGraph* TGraphXclusters;
    TGraph* TGraphYclusters;
    TGraph* TGraphXclusters1;
    TGraph* TGraphYclusters1;
    TGraph* TGraphXclusters2;
    TGraph* TGraphYclusters2;
    TLine anXtrack;
    TLine anYtrack;
    TPaveText* TPaveStat;
    int runId,m_temid;

public:
    EventDisplay(TString="MyRootFile.root", TString geo="",int temid=0);
    virtual ~EventDisplay(){};
    void SetRunId(const int r) { runId = r; }
    void Go(int=-1);


    ClassDef(EventDisplay,1)
};



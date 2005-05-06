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
    ~TriggerEfficiency() {
        delete myEvent;
        delete myTracker;
    }

    void Go(int lastEntry=-1);
};

TriggerEfficiency::TriggerEfficiency(const TString filename,
                                     TString geoFileName) {
    if ( geoFileName.Length() == 0 )
        geoFileName = "$ROOTANALYSISROOT/src/LeaningTower/geometry/TowerBgeometry306000517.txt";
    myTracker = new Tracker;
    // the geometry is only needed to determine the order of the planes
    myTracker->loadGeometry(geoFileName);
    myTracker->SetTower(true);
    myEvent = new Event(filename, myTracker->GetGeometry());
}

void TriggerEfficiency::Go(int lastEntry) {
    int numEntries = myEvent->GetEntries();
    if ( lastEntry == -1 )
        lastEntry = numEntries;
    lastEntry = std::min(numEntries, lastEntry);

    Progress progress;

    UInt_t shouldSum[16];
    UInt_t isSum[16];
    for ( int i=0; i<16; ++i )
        shouldSum[i] = isSum[i] = 0;

#ifdef TRIGGER_EFF_DEBUG
    const TList* myGeometry = myTracker->GetGeometry();
#endif

    for ( int entry=0; entry<lastEntry; ++entry ) { 
        myEvent->Go(entry);
        progress.Go(entry, lastEntry);

        UInt_t digi3row   = myEvent->GetTkrDigi3RowBits();
        UInt_t trgReq3row = myEvent->GetTkrTrgReq3RowBits();
        for ( int i=0; i<16; ++i ) {
	    if ( digi3row&(1<<i) ) // 3-in-a-row hits
                ++shouldSum[i];
            if ( trgReq3row&(1<<i) )
                ++isSum[i];
        }

#ifdef TRIGGER_EFF_DEBUG
        TIter next(myGeometry);
        while ( Layer* plane = (Layer*)next() ) {
            const TString planeName = plane->GetName();
            const int layer = plane->GetLayer();
            if ( (digi3row^trgReq3row)&(1<<layer-1) )
                std::cerr << "eventId " << std::setw(6) << myEvent->GetEventId()
                          << " plane " << std::setw(3) << planeName
                          << " numHits " << myEvent->GetPlaneNumHits(planeName)
                          << " triggerReqOR " << (trgReq3row&(1<<layer-1))
                          << " ToT "<<std::setw(3)<<myEvent->GetToT(planeName,0)
                          << ' ' << std::setw(3) << myEvent->GetToT(planeName,1)
                          << std::endl;
        }
#endif
    }
    
    std::cout << "Layers   shouldTrigger   trigger   efficiency  inefficiency"
              << std::endl;
    for ( int i=0; i<16; ++i ) {
        const float eff = shouldSum[i] > 0 ? 100.*isSum[i]/shouldSum[i] : -100.;
        std::cout << std::setw(2) << i << '-' << std::setw(2) << i+1 << '-'
                  << std::setw(2) << i+2
                  << ' ' << std::setw(13) << shouldSum[i]
                  << ' ' << std::setw(9) << isSum[i]
                  << std::setiosflags(std::ios::fixed)
                  << std::setw(11) << std::setprecision(1) << eff << " % "
                  << std::setw(11) << std::setprecision(1) << 100.-eff << " % "
                  << std::endl;
    }
}

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

TriggerEfficiency::TriggerEfficiency(const TString filename, TString geoFileName) {
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

    const TList* myGeometry = myTracker->GetGeometry();

    Progress progress;

    int shouldTrigger[18];
    int isTrigger[18];
    for ( int l=0; l<18; ++l )
        shouldTrigger[l] = isTrigger[l] = 0;

    for ( int entry=0; entry<lastEntry; ++entry ) { 
        myEvent->Go(entry);
        //        const Int_t eventId = myEvent->GetEventId();
        progress.Go(entry, lastEntry);

        bool hits[18][2];
        bool trigger[18][2];
        for ( int l=0; l<18; ++l )
            for ( int v=0; v<2; ++v )
                hits[l][v] = trigger[l][v] = false;

        TIter next(myGeometry);
        while ( Layer* plane = (Layer*)next() ) {
            const TString planeName = plane->GetName();
            const int TkrNumHits = myEvent->GetPlaneNumHits(planeName);
            const int l = plane->GetLayer();
            const int v = plane->GetView();
            hits[l][v] = TkrNumHits > 0;
            trigger[l][v] = myEvent->GetTriggerReq(planeName,0)
                || myEvent->GetTriggerReq(planeName,1);
            /*
            if ( hits[l][v] ^ trigger[l][v] )
                std::cerr << "eventId " << eventId << " plane " << planeName
                          << ' ' << TkrNumHits << ' ' << trigger[l][v] << std::endl;
            */
        }
        // doing 3-in-a-row
        // index i will measure the efficiency of the layer combination i-1,1,i+1
        for ( int l=1; l<17; ++l )
            if ( hits[l-1][0] && hits[l-1][1] && hits[l][0] && hits[l][1] && hits[l+1][0] && hits[l+1][1] ) { // 3-in-a-row hits
                ++shouldTrigger[l];
                if ( trigger[l-1][0] && trigger[l-1][1] && trigger[l][0] && trigger[l][1] && trigger[l+1][0] && trigger[l+1][1] )
                    ++isTrigger[l];
            }
    }
    std::cout << "Layers   shouldTrigger   trigger   efficiency  inefficiency" << std::endl;
    for ( int l=1; l<17; ++l ) {
        const float eff = shouldTrigger[l] > 0 ? 100. * isTrigger[l] / shouldTrigger[l] : -100.;
        std::cout << std::setw(2) << l-1 << '-' << std::setw(2) << l << '-' << std::setw(2) << l+1
                  << ' ' << std::setw(13) << shouldTrigger[l] << ' ' << std::setw(9) << isTrigger[l]
                  << std::setiosflags(std::ios::fixed)
                  << std::setw(11) << std::setprecision(1) << eff << " % "
                  << std::setw(11) << std::setprecision(1) << 100. - eff << " % "
                  << std::endl;
    }
}

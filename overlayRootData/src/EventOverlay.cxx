//
// The Event class is the top level class for the ASCII to ROOT
// conversion. It contains all the information about a single event.
//

#include "overlayRootData/EventOverlay.h"
#include "Riostream.h"
#include "TCollection.h"  // Declares TIter

ClassImp(EventOverlay)

EventOverlay::EventOverlay() 
{
    m_calOverlayCol  = new TObjArray();
    m_numCalOverlays = -1;

    m_tkrOverlayCol  = new TObjArray();
    m_numTkrOverlays = -1;

    // Assign default values to members
    m_acdOverlayCol  = new TObjArray();
    m_numAcdOverlays = -1;

    Clear();
}

EventOverlay::~EventOverlay() {
  
    m_acdOverlayCol->Delete();
    delete m_acdOverlayCol;
    m_acdOverlayCol = 0;
  
    m_tkrOverlayCol->Delete();
    delete m_tkrOverlayCol;
    m_tkrOverlayCol = 0;
    
    m_calOverlayCol->Delete();
    delete m_calOverlayCol;
    m_calOverlayCol = 0;

    return;
}
  
void EventOverlay::initialize(UInt_t eventId, UInt_t runId, Double_t time, 
                          Double_t liveTime, Bool_t fromMc) 
{
    m_eventId         = eventId;
    m_runId           = runId;
    m_timeStamp       = time;
    m_liveTime        = liveTime;
    m_fromMc          = fromMc;
}

void EventOverlay::Clear(Option_t *option) 
{
    m_eventId = 0;
    m_runId = 0;
    m_timeStamp = 0.0;
    m_liveTime = 0.0;

    m_acdOverlayCol->Clear();
    m_calOverlayCol->Clear();
    m_tkrOverlayCol->Clear();
    m_numAcdOverlays = -1;
    m_numCalOverlays = -1;
    m_numTkrOverlays = -1;
}

void EventOverlay::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Run, Event: " << m_runId << ", " << m_eventId
        << " Time: " << m_timeStamp << " LiveTime: " << m_liveTime << endl;
      
    cout << "Number of CalOverlays " << (m_numCalOverlays+1) << endl;
    if (m_tkrOverlayCol) 
      cout << "Number of TkrOverlays " << m_tkrOverlayCol->GetEntries() << endl;
    else 
      cout << "Number of TkrOverlays 0" << endl;
    cout << "Number of AcdDigis " << m_acdOverlayCol->GetEntries() << endl;
    cout << dec;
}

void EventOverlay::addAcdOverlay(AcdOverlay *acdOverlay)
{
    m_acdOverlayCol->Add(acdOverlay);
    m_numAcdOverlays = m_acdOverlayCol->GetEntries();

    return;
}

const AcdOverlay* EventOverlay::getAcdOverlay(UInt_t i) const 
{
    if (Int_t(i) >= m_acdOverlayCol->GetEntries()) return 0;
    return (AcdOverlay*)m_acdOverlayCol->At(i);
}

void EventOverlay::addTkrOverlay(TkrOverlay *tkrOverlay)
{
    //printf(" add TkrOverlay, %d entries\n",m_tkrOverlayCol->GetEntries());fflush(stdout);
    m_tkrOverlayCol->Add(tkrOverlay);
    m_numTkrOverlays = m_tkrOverlayCol->GetEntries();

    return;
}

const TkrOverlay* EventOverlay::getTkrOverlay(UInt_t i) const 
{
    if (Int_t(i) >= m_tkrOverlayCol->GetEntries()) return 0;
    return (TkrOverlay*)m_tkrOverlayCol->At(i);
}

void EventOverlay::addCalOverlay(CalOverlay* calOverlay) 
{
    // Add the new CalOverlay object to our array
    m_calOverlayCol->Add(calOverlay);
    m_numCalOverlays = m_calOverlayCol->GetEntries();

    return;
}

const CalOverlay* EventOverlay::getCalOverlay(UInt_t i) const 
{
    //if ((m_calOverlayCol) && (m_calOverlayCol->GetEntries()>0)) return (CalOverlay*)m_calOverlayCol->At(i);
    if (Int_t(i) > m_numCalOverlays) return 0;
    return (CalOverlay*)m_calOverlayCol->At(i);
}

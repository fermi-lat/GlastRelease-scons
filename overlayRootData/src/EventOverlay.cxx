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

    m_tkrOverlayCol  = new TObjArray();

    // Assign default values to members
    m_acdOverlayCol  = new TObjArray();

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

    m_acdOverlayCol->Delete();
    m_calOverlayCol->Delete();
    m_tkrOverlayCol->Delete();
    m_acdOverlayCol->Clear();
    m_calOverlayCol->Clear();
    m_tkrOverlayCol->Clear();
}

void EventOverlay::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Run, Event: " << m_runId << ", " << m_eventId
        << " Time: " << m_timeStamp << " LiveTime: " << m_liveTime << endl;
      
    cout << "Number of CalOverlays " << (m_calOverlayCol->GetEntries()) << endl;
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

    return;
}

const AcdOverlay* EventOverlay::getAcdOverlay(UInt_t i) const 
{
    if (Int_t(i) >= m_acdOverlayCol->GetEntries()) return 0;
    return (AcdOverlay*)m_acdOverlayCol->At(i);
}

void EventOverlay::addTkrOverlay(TkrOverlay *tkrOverlay)
{
    m_tkrOverlayCol->Add(tkrOverlay);

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

    return;
}

const CalOverlay* EventOverlay::getCalOverlay(UInt_t i) const 
{
    if (Int_t(i) >= m_calOverlayCol->GetEntries()) return 0;
    return (CalOverlay*)m_calOverlayCol->At(i);
}

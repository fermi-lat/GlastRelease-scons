//
// The Event class is the top level class for the ASCII to ROOT
// conversion. It contains all the information about a single event.
//

#include "digiRootData/DigiEvent.h"
#include <iostream>

ClassImp(DigiEvent)

// Allocate the TClonesArrays and TObjArray just once
TClonesArray *DigiEvent::s_acdDigiStaticCol = 0;
TObjArray *DigiEvent::s_staticTkrDigiCol = 0;
TObjArray *DigiEvent::s_calDigiStaticCol = 0;

DigiEvent::DigiEvent() {

    if (!s_calDigiStaticCol) s_calDigiStaticCol = new TObjArray();
    m_calDigiCol = s_calDigiStaticCol;

    if (!s_staticTkrDigiCol) s_staticTkrDigiCol = new TObjArray();
    m_tkrDigiCol = s_staticTkrDigiCol;

    // Assign default values to members
    if (!s_acdDigiStaticCol) s_acdDigiStaticCol = new TClonesArray("AcdDigi", 1);
    m_acdDigiCol = s_acdDigiStaticCol;
    m_numAcdDigis = -1;

    Clear();
}

DigiEvent::~DigiEvent() {
  
    if(m_acdDigiCol == s_acdDigiStaticCol) s_acdDigiStaticCol = 0;
    m_acdDigiCol->Delete();
    delete m_acdDigiCol;
    m_acdDigiCol = 0;
  
    if(m_tkrDigiCol == s_staticTkrDigiCol) s_staticTkrDigiCol = 0;
    m_tkrDigiCol->Delete();
    delete m_tkrDigiCol;
    m_tkrDigiCol = 0;
    
    if (m_calDigiCol == s_calDigiStaticCol) s_calDigiStaticCol = 0;
    m_calDigiCol->Delete();
    delete m_calDigiCol;
    m_calDigiCol = 0;
}

void DigiEvent::initialize(UInt_t eventId, UInt_t runId, Double_t time, 
                           const L1T& level1, Bool_t fromMc) {
    m_eventId = eventId;
    m_runId = runId;
    m_timeStamp = time;
    m_fromMc = fromMc;
    m_levelOneTrigger = level1;
}

void DigiEvent::Clear(Option_t *option) {
    m_eventId = 0;
    m_runId = 0;
    m_timeStamp = 0.0;
    m_levelOneTrigger.Clear();
    m_calDigiCol->Delete();
    m_tkrDigiCol->Delete();
    m_acdDigiCol->Delete();
    m_numAcdDigis = -1;

}

void DigiEvent::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Run, Event: " << m_runId << ", " << m_eventId
        << " Time: " << m_timeStamp << endl;
    m_levelOneTrigger.Print();
}

AcdDigi* DigiEvent::addAcdDigi(const AcdId& id, Float_t energy, UShort_t *pha, 
                               Bool_t *veto, Bool_t *low, Bool_t *high) {
    // Add a new AcdDigi entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    ++m_numAcdDigis;
    TClonesArray &tiles = *m_acdDigiCol;
    new(tiles[m_numAcdDigis]) AcdDigi(id, energy, pha, veto, low, high);
    return ((AcdDigi*)(tiles[m_numAcdDigis]));
}

const AcdDigi* DigiEvent::getAcdDigi(UShort_t l, UShort_t f, UShort_t r, 
                                     UShort_t c) const {
    // Find a specific AcdDigi in the TClonesArray
    // User supplies a valid AcdDigi identified by layer, face, row, column
    AcdId tempId(l, f, r, c);
    AcdDigi tempDigi = AcdDigi(tempId.getId());
    int index = m_acdDigiCol->BinarySearch(&tempDigi);
    if (index >= 0) return ((AcdDigi*)m_acdDigiCol->At(index));
    return 0;
}

const AcdDigi* DigiEvent::getAcdDigi(const AcdId &id) const {
    // Find a specific AcdDigi in the TClonesArray
    // User supplies a valid AcdId
    AcdDigi tempDigi = AcdDigi(id.getId());
    int index = m_acdDigiCol->BinarySearch(&tempDigi);
    if (index >= 0) return ((AcdDigi*)m_acdDigiCol->At(index));
    return 0;
}

void DigiEvent::addTkrDigi(TkrDigi *digi) {
    m_tkrDigiCol->Add(digi);
}

const TkrDigi* DigiEvent::getTkrDigi(UInt_t i) const {
    if (i >= m_tkrDigiCol->GetEntries()) return 0;
    return (TkrDigi*)m_tkrDigiCol->At(i);
}

void DigiEvent::addCalDigi(CalDigi *cal) {
    m_calDigiCol->Add(cal);
}

const CalDigi* DigiEvent::getCalDigi(UInt_t i) const {
    if (i >= m_calDigiCol->GetEntries()) return 0;
    return (CalDigi*)m_calDigiCol->At(i);
}
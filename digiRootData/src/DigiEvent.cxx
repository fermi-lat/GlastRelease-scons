//
// The Event class is the top level class for the ASCII to ROOT
// conversion. It contains all the information about a single event.
//

#include "digiRootData/DigiEvent.h"
#include <iostream>
#include "TCollection.h"  // Declares TIter

ClassImp(DigiEvent)

// Allocate the TClonesArrays and TObjArray just once
TClonesArray *DigiEvent::s_acdDigiStaticCol = 0;
TObjArray *DigiEvent::s_staticTkrDigiCol = 0;
TClonesArray *DigiEvent::s_calDigiStaticCol = 0;

DigiEvent::DigiEvent() {
    //if (!s_calDigiStaticCol) s_calDigiStaticCol = new TClonesArray("CalDigi",1536);
    if (!s_calDigiStaticCol) s_calDigiStaticCol = new TClonesArray("CalDigi",100);
    m_calDigiCloneCol = s_calDigiStaticCol;
    m_numCalDigis = -1;

    m_calDigiCol = 0;

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
    
    if (m_calDigiCloneCol == s_calDigiStaticCol) s_calDigiStaticCol = 0;
    m_calDigiCloneCol->Delete();
    delete m_calDigiCloneCol;
    m_calDigiCloneCol = 0;

    if (m_calDigiCol) {
        m_calDigiCol->Delete();
        delete m_calDigiCol;
        m_calDigiCol = 0;
    }

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
    const Int_t nd = 10000;
    static Int_t ind=0;
    static TkrDigi* keep[nd];

    m_eventId = 0;
    m_runId = 0;
    m_timeStamp = 0.0;
    m_levelOneTrigger.Clear();
    m_calDigiCloneCol->Clear("C");
    if (m_calDigiCol) {
        m_calDigiCol->Delete();
        delete m_calDigiCol;
        m_calDigiCol = 0;
    }
    //    m_tkrDigiCol->Delete();
    m_acdDigiCol->Clear("C");
    m_numAcdDigis = -1;
    m_numCalDigis = -1;

    //m_tkrDigiCol->Delete();  //<======THIS LINE COMMENTED
   //we delete the objects in keep every nd TkrDigi objects
    //these few lines emulates what TClonesArray is doing

    Int_t n = m_tkrDigiCol->GetEntries();
   for (Int_t i=0;i<n;i++) keep[ind+i] = (TkrDigi*)m_tkrDigiCol->At(i);
   ind += n;
   if (ind > nd-100) {
      for (Int_t j=0;j<ind;j++) delete keep[j];
      ind = 0;
   }
    m_tkrDigiCol->Clear();
}

void DigiEvent::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Run, Event: " << m_runId << ", " << m_eventId
        << " Time: " << m_timeStamp << endl;
    m_levelOneTrigger.Print();
}

AcdDigi* DigiEvent::addAcdDigi(const AcdId& id, const VolumeIdentifier& volId,
                               Float_t energy, UShort_t *pha, 
                               Bool_t *veto, Bool_t *low, Bool_t *high) {
    // Add a new AcdDigi entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    ++m_numAcdDigis;
    TClonesArray &tiles = *m_acdDigiCol;
    new(tiles[m_numAcdDigis]) AcdDigi(id, volId, energy, pha, veto, low, high);
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
    printf(" add TkrDigi, %d entries\n",m_tkrDigiCol->GetEntries());fflush(stdout);
    m_tkrDigiCol->Add(digi);
}

const TkrDigi* DigiEvent::getTkrDigi(UInt_t i) const {
    if (i >= m_tkrDigiCol->GetEntries()) return 0;
    return (TkrDigi*)m_tkrDigiCol->At(i);
}

CalDigi* DigiEvent::addCalDigi() {
    // Add a new CalDigi entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    ++m_numCalDigis;
    TClonesArray &calDigis = *m_calDigiCloneCol;
    new(calDigis[m_numCalDigis]) CalDigi();
    return ((CalDigi*)(calDigis[m_numCalDigis]));
}

const CalDigi* DigiEvent::getCalDigi(UInt_t i) const {
    if (m_calDigiCol) return (CalDigi*)m_calDigiCol->At(i);
    if (i > m_numCalDigis) return 0;
    return (CalDigi*)m_calDigiCloneCol->At(i);
}

const TClonesArray* DigiEvent::getCalDigiCol() { 
    // Purpose and Method:  Provide access to the collection of CalDigis.
    // For backward compatibilty, we must handle the case where the m_calDigiCol
    // (a TObjArray*) exists - in that case, we fill up the TClonesArray to allow the
    // user to reap the benefits of the TClonesArray, and to provide a reasonable 
    // method for calling this routine and receiving a non-null return value.
    // This conversion to a TClonesArray will only happen the first time this method is
    // called for a particular event.

    // Maintain Backward compatibility
    if ((m_calDigiCol) && (m_calDigiCloneCol->GetEntries()==0)) {
        TIter objIter(m_calDigiCol);
        CalDigi *calDigiObj = 0;
        while (calDigiObj = (CalDigi*)objIter.Next()) {
            CalDigi *newObj = addCalDigi();
            newObj->initialize(calDigiObj->getMode(), calDigiObj->getPackedId());
            if (calDigiObj->getMode() == CalXtalId::BESTRANGE) {
                const CalXtalReadout *readout = calDigiObj->getXtalReadout(0);
                Char_t rangePlus = readout->getRange(CalXtalId::POS);
                UInt_t adcPlus = readout->getAdc(CalXtalId::POS);
                Char_t rangeMin = readout->getRange(CalXtalId::NEG);
                UInt_t adcMin = readout->getAdc(CalXtalId::NEG);
                newObj->addReadout(rangePlus, adcPlus, rangeMin, adcMin);
            } else { // Handle AllRange
                Int_t range;
                for (range = CalXtalId::LEX8; range <= CalXtalId::HEX1; range++) {
                    const CalXtalReadout *readout = calDigiObj->getXtalReadout(range);
                    Char_t rangePlus = readout->getRange(CalXtalId::POS);
                    UInt_t adcPlus = readout->getAdc(CalXtalId::POS);
                    Char_t rangeMin = readout->getRange(CalXtalId::NEG);
                    UInt_t adcMin = readout->getAdc(CalXtalId::NEG);
                    newObj->addReadout(rangePlus, adcPlus, rangeMin, adcMin);
                }
            }
        }
    }
    // Now return the TClonesArray*
    return m_calDigiCloneCol; 
}

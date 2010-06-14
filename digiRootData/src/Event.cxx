///////////////////////////////////////////////////////////////////////////
//
// The Event class is the top level class for the ASCII to ROOT
// conversion. It contains all the information about a single event.
// It contains TObjArrays (typically one per detector component) of
// information about hit detector elements, as well as functions to
// access these arrays.
//
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/ACDTile.h"    // So that we can manually delete ACDTiles
                              // [ See Event::Clean() for details ]
#include "TClass.h"
#include "digiRootData/Event.h"
#include "TMath.h"

ClassImp(Event)
//__________________________________________________________________________
Event::Event() {
    // Default constructor.
    // Assign default values to members
    m_ACD = 0;
    m_CAL = 0;
    m_TKR = 0;
    m_MCPart = 0;
    m_L1Trigger = 0;
    m_Tagger=0;  //THB???
    //    m_HSK = 0;
}
//_________________________________________________________________________
Event::~Event() {
  // Destructor
  Clean();
}
//__________________________________________________________________________
void Event::Create() {
    m_ACD = new TObjArray();
    m_CAL = new TObjArray();
    m_TKR = new TObjArray();
    m_MCPart = new TObjArray();
    m_L1Trigger = 0;
    //    m_HSK = 0;
}
//_________________________________________________________________________
void Event::Clean() {
    // Clears and destroys all TObjArrays contained in the Event object.
    // Necessary both for destruction and between subsequent fills and
    // writes/reads to/from file.

    if (m_CAL) {
        int nEntries = m_CAL->GetEntries();
        for (int i=0; i<nEntries; i++)
            delete m_CAL->At(i);
        m_CAL->Clear();
        delete m_CAL;
        m_CAL = 0;
    }

    if (m_ACD) {
        int nEntries = m_ACD->GetEntries();
        for (int i=0; i<nEntries; i++)
            delete m_ACD->At(i);
        m_ACD->Clear();
        delete m_ACD;
        m_ACD = 0;
    }

    if (m_TKR) {
        int nEntries = m_TKR->GetEntries();
        for (int i=0; i<nEntries; i++)
            delete m_TKR->At(i);
        m_TKR->Clear();
        delete m_TKR;
        m_TKR = 0;
    }

    if (m_MCPart) {
        int nEntries = m_MCPart->GetEntries();
        for (int i=0; i<nEntries; i++)
            delete m_MCPart->At(i);
        m_MCPart->Clear();
        delete m_MCPart;
        m_MCPart = 0;
    }

    if (m_L1Trigger) {
        delete m_L1Trigger;
        m_L1Trigger = 0;
    }
    /*
    if (m_HSK) {
      delete m_HSK;
      m_HSK = 0;
    }
    */
    return;
}

//________________________________________________________________________
TObjArray* Event::getACD() {
    // Provides access to list of ACD tiles
    return m_ACD;
}
//________________________________________________________________________
TObjArray* Event::getCAL(){
    // Provides access to list of Calorimeter CsI logs
    return m_CAL;
}
//________________________________________________________________________
TObjArray* Event::getMCPart(){
    // Provides access to list of MC particles from
    // reconstruction of this event
    return m_MCPart;
}
//________________________________________________________________________
TObjArray* Event::getTKR(){
    // Provides access to list of tracker layers
    return m_TKR;
}
//________________________________________________________________________
/// Implement the streamer ourselves for now...to take advantage
/// of schema evolution - and to allow our TBEvent class to handle
/// both old (<= Root v2.25) Root files, and new (>= Root v3.00) files
void Event::Streamer(TBuffer &R__b)
{
   // Stream an object of class Event.
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); 
      if (R__v > 1) 
      { 
          Event::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
          return;
      }
      /// Old Versions
      TObject::Streamer(R__b);
      R__b >> m_run;
      R__b >> m_event;
      R__b >> m_CAL;
      R__b >> m_ACD;
      R__b >> m_TKR;
      R__b >> m_MCPart;
      R__b >> m_L1Trigger;
      R__b >> m_Tagger;
   } else {
       Event::Class()->WriteBuffer(R__b, this);
   }
}
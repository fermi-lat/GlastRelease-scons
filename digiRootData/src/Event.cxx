///////////////////////////////////////////////////////////////////////////
//
// The Event class is the top level class for the ASCII to ROOT
// conversion. It contains all the information about a single event.
//
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/Event.h"

ClassImp(Event)

// Allocate the TClonesArrays and TObjArray just once
TClonesArray *Event::m_staticAcdDigiVec = 0;
TClonesArray *Event::m_staticCalDigiVec = 0;
TObjArray *Event::m_staticTkrDigiVec = 0;

//__________________________________________________________________________
Event::Event() {
    // Default constructor.
    // Assign default values to members
    if (!m_staticAcdDigiVec) m_staticAcdDigiVec = new TClonesArray("AcdDigi", 24);
    m_AcdDigiVec = m_staticAcdDigiVec;
    m_numTiles = -1;

    if (!m_staticCalDigiVec) m_staticCalDigiVec = new TClonesArray("CalDigi", 80);
    m_CalDigiVec = m_staticCalDigiVec;
    m_numLogs = -1;

    if (!m_staticTkrDigiVec) m_staticTkrDigiVec = new TObjArray();
    m_TkrDigiVec = m_staticTkrDigiVec;
    m_numLayers = -1;

        
  //  m_Tagger = 0;  
    m_eventId = 0;
    m_run = 0;
}
//_________________________________________________________________________
Event::~Event() {
  // Destructor    
    if(m_AcdDigiVec == m_staticAcdDigiVec) m_staticAcdDigiVec = 0;
    m_AcdDigiVec->Delete();
    delete m_AcdDigiVec;
    m_AcdDigiVec = 0;
  
    if(m_CalDigiVec == m_staticCalDigiVec) m_staticCalDigiVec = 0;
    m_CalDigiVec->Delete();
    delete m_CalDigiVec;
    m_CalDigiVec = 0;

    if(m_TkrDigiVec == m_staticTkrDigiVec) m_staticTkrDigiVec = 0;
    m_TkrDigiVec->Delete();
    delete m_TkrDigiVec;
    m_TkrDigiVec = 0;
    
}
//_________________________________________________________________________
void Event::Clean(Option_t *option) {
    m_AcdHeader.Clean();
    m_AcdDigiVec->Clear();
    m_numTiles = -1;

    m_CalHeader.Clean();
    m_CalDigiVec->Clear();
    m_numLogs = -1;

    m_TkrHeader.Clean();
    m_TkrDigiVec->Delete();
    m_numLayers = -1;

    m_L1T.Clean();
}


AcdDigi* Event::addAcdDigi(UShort_t id) {
    // Add a new AcdDigi entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    ++m_numTiles;
    TClonesArray &tiles = *m_AcdDigiVec;
    new(tiles[m_numTiles]) AcdDigi(id);
    return ((AcdDigi*)(tiles[m_numTiles]));
}

CalDigi* Event::addCalDigi() {
    // Add a new CalDigi entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    ++m_numLogs;
    TClonesArray &logs = *m_CalDigiVec;
    new(logs[m_numLogs]) CalDigi();
    return ((CalDigi*)(logs[m_numLogs]));
}

void Event::addTkrDigi(TkrDigi *layer) {
    m_TkrDigiVec->Add(layer);
    ++m_numLayers;
}

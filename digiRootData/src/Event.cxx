///////////////////////////////////////////////////////////////////////////
//
// The Event class is the top level class for the ASCII to ROOT
// conversion. It contains all the information about a single event.
//
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/Event.h"

ClassImp(Event)

// Allocate the TClonesArrays and TObjArray just once
TClonesArray *Event::m_staticAcdData = 0;
TClonesArray *Event::m_staticCalData = 0;
TObjArray *Event::m_staticTkrData = 0;

//__________________________________________________________________________
Event::Event() {
    // Default constructor.
    // Assign default values to members
    if (!m_staticAcdData) m_staticAcdData = new TClonesArray("ACDTile", 24);
    m_AcdData = m_staticAcdData;
    m_numTiles = -1;

    if (!m_staticCalData) m_staticCalData = new TClonesArray("CalHit", 80);
    m_CalData = m_staticCalData;
    m_numLogs = -1;

    if (!m_staticTkrData) m_staticTkrData = new TObjArray();
    m_TkrData = m_staticTkrData;
    m_numLayers = -1;

        
  //  m_Tagger = 0;  
    m_eventId = 0;
    m_run = 0;
}
//_________________________________________________________________________
Event::~Event() {
  // Destructor    
    if(m_AcdData == m_staticAcdData) m_staticAcdData = 0;
    m_AcdData->Delete();
    delete m_AcdData;
    m_AcdData = 0;
  
    if(m_CalData == m_staticCalData) m_staticCalData = 0;
    m_CalData->Delete();
    delete m_CalData;
    m_CalData = 0;

    if(m_TkrData == m_staticTkrData) m_staticTkrData = 0;
    m_TkrData->Delete();
    delete m_TkrData;
    m_TkrData = 0;
    
}
//_________________________________________________________________________
void Event::Clean(Option_t *option) {
    m_AcdHeader.Clean();
    m_AcdData->Clear();
    m_numTiles = -1;

    m_CalHeader.Clean();
    m_CalData->Clear();
    m_numLogs = -1;

    m_TkrHeader.Clean();
    m_TkrData->Delete();
    m_numLayers = -1;

    m_L1T.Clean();
}


ACDTile* Event::AddTile(UShort_t id) {
    // Add a new ACDTile entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    ++m_numTiles;
    TClonesArray &tiles = *m_AcdData;
    new(tiles[m_numTiles]) ACDTile(id);
    return ((ACDTile*)(tiles[m_numTiles]));
}

CalHit* Event::AddLog() {
    // Add a new ACDTile entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    ++m_numLogs;
    TClonesArray &logs = *m_CalData;
    new(logs[m_numLogs]) CalHit();
    return ((CalHit*)(logs[m_numLogs]));
}

void Event::AddLayer(TkrLayer *layer) {
    m_TkrData->Add(layer);
    ++m_numLayers;
}

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
TClonesArray *Event::m_staticXgtDigiVec = 0;
TClonesArray *Event::m_staticCalDigiVec = 0;
TObjArray *Event::m_staticTkrDigiVec = 0;

//__________________________________________________________________________
Event::Event() {
    // Default constructor.
    // Assign default values to members
    if (!m_staticAcdDigiVec) m_staticAcdDigiVec = new TClonesArray("AcdTile", 24);
    m_AcdDigiVec = m_staticAcdDigiVec;
    m_numTiles = -1;

    if (!m_staticXgtDigiVec) m_staticXgtDigiVec = new TClonesArray("AcdTile", 4);
    m_XgtDigiVec = m_staticXgtDigiVec;
    m_numXgts = -1;

    if (!m_staticCalDigiVec) m_staticCalDigiVec = new TClonesArray("CalLog", 80);
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

    if(m_XgtDigiVec == m_staticXgtDigiVec) m_staticXgtDigiVec = 0;
    m_XgtDigiVec->Delete();
    delete m_XgtDigiVec;
    m_XgtDigiVec = 0;
  
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

    m_XgtDigiVec->Clear();
    m_numXgts = -1;

    m_CalHeader.Clean();
    m_CalDigiVec->Clear();
    m_numLogs = -1;

    m_TkrHeader.Clean();
    m_TkrDigiVec->Delete();
    m_numLayers = -1;

    m_L1T.Clean();
}


AcdTile* Event::addAcdTile(Short_t id) {
    // Add a new AcdTile entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    ++m_numTiles;
    TClonesArray &tiles = *m_AcdDigiVec;
    new(tiles[m_numTiles]) AcdTile(id);
    return ((AcdTile*)(tiles[m_numTiles]));
}

const AcdTile* Event::getAcdTile(UShort_t id) {
    // Find a specific AcdTile in the TClonesArray
    // User supplies a valid AcdTile id
    AcdTile tempTile = AcdTile(id);
    int index = m_AcdDigiVec->BinarySearch(&tempTile);
    if (index >= 0) return ((AcdTile*)m_AcdDigiVec->At(index));
    return 0;
}

AcdTile* Event::addXgt(Short_t id) {
    // Add a new XGT(AcdTile) entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    ++m_numXgts;
    TClonesArray &xgts = *m_XgtDigiVec;
    new(xgts[m_numXgts]) AcdTile(id);
    return ((AcdTile*)(xgts[m_numXgts]));
}

const AcdTile* Event::getXgt(UShort_t id) {
    // Find a specific XGT(AcdTile) in the TClonesArray
    // User supplies a valid AcdTile id
    AcdTile tempTile = AcdTile(id);
    int index = m_XgtDigiVec->BinarySearch(&tempTile);
    if (index >= 0) return ((AcdTile*)m_XgtDigiVec->At(index));
    return 0;
}


CalLog* Event::addCalLog() {
    // Add a new CalLog entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    ++m_numLogs;
    TClonesArray &logs = *m_CalDigiVec;
    new(logs[m_numLogs]) CalLog();
    return ((CalLog*)(logs[m_numLogs]));
}

void Event::addTkrLayer(TkrLayer *layer) {
    m_TkrDigiVec->Add(layer);
    ++m_numLayers;
}

const TkrLayer* Event::getTkrLayer(unsigned int layerNum) {
    return ((TkrLayer*)m_TkrDigiVec->At(layerNum));
}

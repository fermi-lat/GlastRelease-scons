//
// The Event class is the top level class for the ASCII to ROOT
// conversion. It contains all the information about a single event.
//

#include "digiRootData/DigiEvent.h"

ClassImp(DigiEvent)

// Allocate the TClonesArrays and TObjArray just once
//TClonesArray *DigiEvent::m_staticAcdDigiVec = 0;
TObjArray *DigiEvent::s_staticTkrDigiCol = 0;
TObjArray *DigiEvent::s_calDigiStaticCol = 0;

DigiEvent::DigiEvent() {

    if (!s_calDigiStaticCol) s_calDigiStaticCol = new TObjArray();
    m_calDigiCol = s_calDigiStaticCol;

    if (!s_staticTkrDigiCol) s_staticTkrDigiCol = new TObjArray();
    m_tkrDigiCol = s_staticTkrDigiCol;

    Clear();

    /*
    // Assign default values to members
    if (!m_staticAcdDigiVec) m_staticAcdDigiVec = new TClonesArray("AcdTile", 24);
    m_AcdDigiVec = m_staticAcdDigiVec;
    m_numTiles = -1;

    if (!m_staticTkrDigiVec) m_staticTkrDigiVec = new TObjArray();
    m_TkrDigiVec = m_staticTkrDigiVec;
    m_numLayers = -1;
*/
}

DigiEvent::~DigiEvent() {
  /*
    if(m_AcdDigiVec == m_staticAcdDigiVec) m_staticAcdDigiVec = 0;
    m_AcdDigiVec->Delete();
    delete m_AcdDigiVec;
    m_AcdDigiVec = 0;
  */
    if(m_tkrDigiCol == s_staticTkrDigiCol) s_staticTkrDigiCol = 0;
    m_tkrDigiCol->Delete();
    delete m_tkrDigiCol;
    m_tkrDigiCol = 0;
    
    if (m_calDigiCol == s_calDigiStaticCol) s_calDigiStaticCol = 0;
    m_calDigiCol->Delete();
    delete m_calDigiCol;
    m_calDigiCol = 0;
}

void DigiEvent::initialize(UInt_t eventId, UInt_t runId, Bool_t fromMc) {
    m_eventId = eventId;
    m_runId = runId;
    m_fromMc = fromMc;
}

void DigiEvent::Clear(Option_t *option) {
    m_eventId = 0;
    m_runId = 0;
    m_calDigiCol->Delete();
    m_tkrDigiCol->Delete();

    /*
    m_AcdHeader.Clean();
    m_AcdDigiVec->Clear();
    m_numTiles = -1;

    m_CalHeader.Clean();

    m_TkrHeader.Clean();

    m_L1T.Clean();

    m_liveTime.Clean();
    */
}

void DigiEvent::Print(Option_t *option) const {

}

/*
AcdTile* DigiEvent::addAcdTile(UInt_t id, short base, short used) {
    // Add a new AcdTile entry, note that
    // TClonesArrays can only be filled via
    // a new with placement call
    ++m_numTiles;
    TClonesArray &tiles = *m_AcdDigiVec;
    new(tiles[m_numTiles]) AcdTile(id, base, used);
    return ((AcdTile*)(tiles[m_numTiles]));
}

const AcdTile* DigiEvent::getAcdTile(UInt_t id) {
    // Find a specific AcdTile in the TClonesArray
    // User supplies a valid AcdTile id, in base 10 or base 2
    // default is base 10
    AcdTile tempTile = AcdTile(id);
    int index = m_AcdDigiVec->BinarySearch(&tempTile);
    if (index >= 0) return ((AcdTile*)m_AcdDigiVec->At(index));
    return 0;
}

const AcdTile* DigiEvent::getAcdTile(short l, short f, short r, short c) {
    // Find a specific AcdTile in the TClonesArray
    // User supplies a valid AcdTile identified by layer, face, row, column
    AcdId tempId(l, f, r, c);
    AcdTile tempTile = AcdTile(tempId.getId());
    int index = m_AcdDigiVec->BinarySearch(&tempTile);
    if (index >= 0) return ((AcdTile*)m_AcdDigiVec->At(index));
    return 0;
}

const AcdTile* DigiEvent::getAcdTile(AcdId &id) {
    // Find a specific AcdTile in the TClonesArray
    // User supplies a valid AcdId
    AcdTile tempTile = AcdTile(id.getId());
    int index = m_AcdDigiVec->BinarySearch(&tempTile);
    if (index >= 0) return ((AcdTile*)m_AcdDigiVec->At(index));
    return 0;
}



const TkrLayer* DigiEvent::getTkrLayer(unsigned int layerNum) {
    TkrLayer tempLayer = TkrLayer();
    unsigned int planeNum = layerNum >> 1;
    tempLayer.setPlaneNum(planeNum);
    if ( (planeNum % 2) == 0 ) {    // even plane
        tempLayer.setXY((layerNum&1) ? TkrLayer::X : TkrLayer::Y);
    } else {                        // odd plane
        tempLayer.setXY((layerNum&1) ? TkrLayer::Y : TkrLayer::X);
    }
    int index = m_TkrDigiVec->BinarySearch(&tempLayer);
    if (index >= 0) return ((TkrLayer*)m_TkrDigiVec->At(index));
    return 0;
}
*/

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
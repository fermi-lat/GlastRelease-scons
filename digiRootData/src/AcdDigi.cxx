#include "digiRootData/AcdDigi.h"


ClassImp(AcdDigi)

TClonesArray *AcdDigi::m_staticAcdArray = 0;


AcdDigi::AcdDigi() {
    if (!m_staticAcdArray) m_staticAcdArray = new TClonesArray("AcdTile", 24);
    m_tiles = m_staticAcdArray;
    m_numTiles = -1;
    m_eventId = 0;
    m_timerWord = 0;
    m_TREQ_VETO_status = 0;
    m_deadTime = 0;
    m_deadTimeCause = 0;
}

AcdDigi::~AcdDigi() {
    // Destructor
    if (m_tiles == m_staticAcdArray) m_staticAcdArray = 0;
    m_tiles->Delete();
    delete m_tiles;
}


AcdTile* AcdDigi::Add(TileID *id) {
    ++m_numTiles;    
    TClonesArray &tiles = *m_tiles;
    new(tiles[m_numTiles]) AcdTile(id);
    return ((AcdTile*)(tiles[m_numTiles]));

}

void AcdDigi::Clean(Option_t *option) {
    // Call Delete rather than Clear
    // since the elements in this array
    // have allocated memory
    m_tiles->Delete(option);
}
///////////////////////////////////////////////////////////////////////////
//									 
// The ACDTile class contains the information about a single ACD tile.
// This includes the PHA value, and above thresh information.  
//									 
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/ACDTile.h"

ClassImp(ACDTile)

///________________________________________________________________________
ACDTile::ACDTile() : m_tag(0){
    // Default constructor
    m_tileID.setID(0);
}
///________________________________________________________________________
//ACDTile::ACDTile(TileID& id) : m_tag(0) {
    // Create a ACDTile object with tileID of id
//    m_tileID = id;
//}
ACDTile::ACDTile(UShort_t i) : m_tag(0){//, m_tileId(i) {
    m_tileID.setID(i);
}
//_________________________________________________________________________
ACDTile::~ACDTile(){
    // Destructor 
   // if (m_tileID) {
   //	  delete m_tileID;
   //	  m_tileID = 0;
    //}
}
//_________________________________________________________________________
Int_t ACDTile::Compare(const TObject *obj) const {
    if (this == obj) return 0;
    if (ACDTile::Class() != obj->IsA()) return -1;
    UInt_t id_this = m_tileID.getID();
    UInt_t id_tile = ((ACDTile*)obj)->getID()->getID();
    if (id_this == id_tile)
	return 0;
    else
	return (id_this > id_tile) ? 1 : -1; 
}
//_________________________________________________________________________
Bool_t ACDTile::IsSortable() const {
  return kTRUE;
}
//_________________________________________________________________________
UChar_t ACDTile::getHit() {
    // Returns 1 if tile above thresh, 0 otherwise
    return ((m_tag >> ACD_V_HIT) & ACD_M_HIT);
}
//_________________________________________________________________________
UShort_t ACDTile::getPHA() {
    // Returns the PMT value for this tile
    return ((m_tag >> ACD_V_PMT) & ACD_M_PMT);
}
//_________________________________________________________________________
Bool_t ACDTile::setPHA(UShort_t pmtVal)
{
    // Set the PMT/PHA value for this tile
    // Valid PMT/PHA values are in the range of [0,2047]
    if (pmtVal & ~ACD_M_PMT)
	return kFALSE;
    else {
	m_tag &= ~(ACD_M_PMT << ACD_V_PMT);
	m_tag |= pmtVal;
	return kTRUE;
    }
}
//_________________________________________________________________________
Bool_t ACDTile::setHit(UChar_t hitVal)
{
    // Set the 'hit' status of this tile
    // 1 == above thresh, 0 == below
    if (hitVal & ~ACD_M_HIT) 
	return kFALSE;
    else {
	m_tag &= ~(ACD_M_HIT << ACD_V_HIT);
	m_tag |= hitVal;
	return kTRUE;
    }
}

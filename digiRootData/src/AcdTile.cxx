///////////////////////////////////////////////////////////////////////////
//									 
// The AcdTile class contains the information about a single ACD tile.
// This includes the PHA value, and above thresh information.  
//									 
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/AcdTile.h"

ClassImp(AcdTile)

///________________________________________________________________________
AcdTile::AcdTile() : m_tag(0){
    // Default constructor
    m_tileId.setId(0);
}
///________________________________________________________________________
AcdTile::AcdTile(UShort_t i) : m_tag(0){//, m_tileId(i) {
    m_tileId.setId(i);
}
//_________________________________________________________________________
AcdTile::~AcdTile(){
    // Destructor 
}
//_________________________________________________________________________
Int_t AcdTile::Compare(const TObject *obj) const {
    if (this == obj) return 0;
    if (AcdTile::Class() != obj->IsA()) return -1;
    UInt_t id_this = m_tileId.getId();
    UInt_t id_tile = ((AcdTile*)obj)->getId()->getId();
    if (id_this == id_tile)
	return 0;
    else
	return (id_this > id_tile) ? 1 : -1; 
}
//_________________________________________________________________________
Bool_t AcdTile::IsSortable() const {
  return kTRUE;
}
//_________________________________________________________________________
UChar_t AcdTile::getVeto() {
    // Returns 1 if tile above thresh, 0 otherwise
    return ((m_tag >> ACD_V_HIT) & ACD_M_HIT);
}
//_________________________________________________________________________
UShort_t AcdTile::getPulseHeight() {
    // Returns the PMT value for this tile
    return ((m_tag >> ACD_V_PMT) & ACD_M_PMT);
}
//_________________________________________________________________________
Bool_t AcdTile::setPulseHeight(UShort_t pmtVal)
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
Bool_t AcdTile::setVeto(UChar_t hitVal)
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

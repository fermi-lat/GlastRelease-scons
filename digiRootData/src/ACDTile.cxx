///////////////////////////////////////////////////////////////////////////
//                                                                       
// The ACDTile class contains the information about a single ACD tile.
// This includes the PHA value, and above thresh information.  
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/ACDTile.h"
#include "TClass.h"

ClassImp(ACDTile)

///________________________________________________________________________
ACDTile::ACDTile() : m_tag(0){
    // Default constructor
    m_tileID = 0;
}
///________________________________________________________________________
ACDTile::ACDTile(TileID* id) : m_tag(0) {
    // Create a ACDTile object with tileID of id
    m_tileID = id;
}
//_________________________________________________________________________
ACDTile::~ACDTile(){
    // Destructor 
    if (m_tileID)
        delete m_tileID;
}
//_________________________________________________________________________
Int_t ACDTile::Compare(const TObject *obj) const {
    if (this == obj) return 0;
    if (ACDTile::Class() != obj->IsA()) return -1;
    UInt_t id_this = m_tileID->getID();
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
UShort_t ACDTile::getPMT() {
    // Returns the PMT value for this tile
    return ((m_tag >> ACD_V_PMT) & ACD_M_PMT);
}
//_________________________________________________________________________
Bool_t ACDTile::setPMT(UShort_t pmtVal)
{
    // Set the PMT value for this tile
    // Valid PMT values are in the range of [0,2047]
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
//_________________________________________________________________________
/// Implement the streamer ourselves for now...to take advantage
/// of schema evolution - and to allow our TBEvent class to handle
/// both old (<= Root v2.25) Root files, and new (>= Root v3.00) files
void ACDTile::Streamer(TBuffer &R__b)
{
   // Stream an object of class ACDTile.
   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); 
      if (R__v > 1) { 
          ACDTile::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
          return;
      }
      TObject::Streamer(R__b);
      R__b >> m_tag;
      R__b >> m_tileID;
   } else {
       ACDTile::Class()->WriteBuffer(R__b, this);
   }
}
//									 
// The AcdDigi class contains the information about a single ACD tile.
// This includes the PHA value, and above thresh information.  
//									 

#include "digiRootData/AcdDigi.h"
#include <iostream>

ClassImp(AcdDigi)

AcdDigi::AcdDigi() {
    m_packed[AcdDigi::A] = 0;
    m_packed[AcdDigi::B] = 0;
    m_id.setId(0, 1, 2);
}

AcdDigi::AcdDigi(const AcdId& id) : m_id(id) {
    m_packed[AcdDigi::A] = 0;
    m_packed[AcdDigi::B] = 0;
}

AcdDigi::AcdDigi(const AcdId& id, Float_t energy, UShort_t *pha,
                 Bool_t *veto, Bool_t *low, Bool_t *high) 
                 : m_id(id)
{
    AcdDigi::PmtId pmt = AcdDigi::A;
    initPackedWord(pmt, pha[pmt], veto[pmt], low[pmt], high[pmt]); 
    pmt = AcdDigi::B;
    initPackedWord(pmt, pha[pmt], veto[pmt], low[pmt], high[pmt]); 

    m_energy = energy;
}


void AcdDigi::Clear(Option_t *option) {

}

void AcdDigi::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Id: " << endl;
    m_id.Print();
    cout << "PMT A: " << getPulseHeight(AcdDigi::A)  
        << " (low, veto, high): (" << Short_t(getLowDiscrim(AcdDigi::A)) << ", "
        << Short_t(getVeto(AcdDigi::A)) << ", " 
        << Short_t(getHighDiscrim(AcdDigi::A)) << ")" << endl;
    cout << "PMT B: " << getPulseHeight(AcdDigi::B) 
        << " (low, veto, high): (" << Short_t(getLowDiscrim(AcdDigi::B)) << ", "
        << Short_t(getVeto(AcdDigi::B)) << ", " 
        << Short_t(getHighDiscrim(AcdDigi::B)) << ")" << endl;
}


Bool_t AcdDigi::getLowDiscrim(AcdDigi::PmtId pmt) const {
    return ( (m_packed[pmt] >> LOW_SHIFT) & 1 );
}

Bool_t AcdDigi::getHighDiscrim(AcdDigi::PmtId pmt) const {

    return ( (m_packed[pmt] >> HIGH_SHIFT) & 1 );
}

Bool_t AcdDigi::getVeto(AcdDigi::PmtId pmt) const {
    return ( (m_packed[pmt] >> VETO_SHIFT) & 1 );
}

UShort_t AcdDigi::getPulseHeight(AcdDigi::PmtId pmt) const {
    return ((m_packed[pmt] >> PMT_SHIFT) & PMT_MAX);
}

void AcdDigi::initPackedWord(AcdDigi::PmtId pmt, UShort_t pha, Bool_t veto,
                    Bool_t low, Bool_t high) 
{
    m_packed[pmt] = 0;
    if (high == kTRUE) m_packed[pmt] |= (1 << HIGH_SHIFT);
    if (veto == kTRUE) m_packed[pmt] |= (1 << VETO_SHIFT);
    if (low == kTRUE) m_packed[pmt] |= (1 << LOW_SHIFT);
    
    if (pha & ~PMT_MAX) {
    } else {
        m_packed[pmt] |= pha;
    }
}

Int_t AcdDigi::Compare(const TObject *obj) const {
    const short base = 2;
    if (this == obj) return 0;
    if (AcdDigi::Class() != obj->IsA()) return -1;
    UInt_t id_this = m_id.getId(base);
    UInt_t id_tile = ((AcdDigi*)obj)->getId().getId(base);
    if (id_this == id_tile)
	return 0;
    else
	return (id_this > id_tile) ? 1 : -1; 
}

Bool_t AcdDigi::IsSortable() const {
  return kTRUE;
}


//									 
// The AcdDigi class contains the information about a single ACD tile.
// This includes the PHA value, and above thresh information.  
//									 

#include "digiRootData/AcdDigi.h"
#include "Riostream.h"

ClassImp(AcdDigi)

AcdDigi::AcdDigi() {
    m_packed[AcdDigi::A] = 0;
    m_packed[AcdDigi::B] = 0;
    m_id.setId(0, 1, 2);
    m_packedLdf[AcdDigi::A] = 0;
    m_packedLdf[AcdDigi::B] = 0;
    m_tileName = TString("");
    m_tileNumber = -1;
}

AcdDigi::AcdDigi(const AcdId& id) : m_id(id) {
    m_packed[AcdDigi::A] = 0;
    m_packed[AcdDigi::B] = 0;
    m_packedLdf[AcdDigi::A] = 0;
    m_packedLdf[AcdDigi::B] = 0;
    m_tileName = TString("");
    m_tileNumber = -1;
}

AcdDigi::AcdDigi(const AcdId& id, const VolumeIdentifier& volId, Float_t energy, 
                 UShort_t *pha, Bool_t *veto, Bool_t *low, Bool_t *high) 
                 : m_id(id), m_volId(volId), m_tileNumber(-1), m_tileName("")
{
    AcdDigi::PmtId pmt = AcdDigi::A;
    initPackedWord(pmt, pha[pmt], veto[pmt], low[pmt], high[pmt]); 
    initPackedLdfWord(pmt, LOW, NOERROR, NOERROR);
    pmt = AcdDigi::B;
    initPackedWord(pmt, pha[pmt], veto[pmt], low[pmt], high[pmt]); 
    initPackedLdfWord(pmt, LOW, NOERROR, NOERROR);

    m_energy = energy;
}

void AcdDigi::initLdfParameters(const char *name, Int_t number, Range *rangeCol, ParityError *errCol, ParityError *headerParityCol) {
    m_tileName = name;
    m_tileNumber = number;
    AcdDigi::PmtId pmt = AcdDigi::A;
    initPackedLdfWord(pmt, rangeCol[pmt], errCol[pmt], headerParityCol[pmt]);
    pmt = AcdDigi::B;
    initPackedLdfWord(pmt, rangeCol[pmt], errCol[pmt], headerParityCol[pmt]);
}

void AcdDigi::Clear(Option_t *option) {

}

void AcdDigi::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout.precision(2);
    cout << "Id: " << endl;
    m_id.Print();
    cout << "name, number: " << m_tileName.Data() << ", " << m_tileNumber << endl;
    cout << "PMT A: " << getPulseHeight(AcdDigi::A)  
        << " (low, veto, high): (" << Short_t(getLowDiscrim(AcdDigi::A)) << ", "
        << Short_t(getVeto(AcdDigi::A)) << ", " 
        << Short_t(getHighDiscrim(AcdDigi::A)) << ")" << endl;
    cout << " (range, error): (" << Short_t(getRange(AcdDigi::A)) << ", "
        << Short_t(getOddParityError(AcdDigi::A)) << ")" << endl;
    cout << "PMT B: " << getPulseHeight(AcdDigi::B) 
        << " (low, veto, high): (" << Short_t(getLowDiscrim(AcdDigi::B)) << ", "
        << Short_t(getVeto(AcdDigi::B)) << ", " 
        << Short_t(getHighDiscrim(AcdDigi::B)) << ")" << endl;
    cout  << " (range, error): (" << Short_t(getRange(AcdDigi::B)) << ", "
        << Short_t(getOddParityError(AcdDigi::B)) << ")" << endl;
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

AcdDigi::Range AcdDigi::getRange(AcdDigi::PmtId pmt) const {
    return (((m_packedLdf[pmt] >> RANGE_SHIFT) & 1 ) ? HIGH : LOW);
}

AcdDigi::ParityError AcdDigi::getOddParityError(AcdDigi::PmtId pmt) const {
    return (((m_packedLdf[pmt] >> ERROR_SHIFT) & 1 ) ? ERROR : NOERROR);
}

AcdDigi::ParityError AcdDigi::getHeaderParityError(AcdDigi::PmtId pmt) const {
    return (((m_packedLdf[pmt] >> HEADERPARITY_SHIFT) & 1 ) ? ERROR : NOERROR);
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

void AcdDigi::initPackedLdfWord(AcdDigi::PmtId pmt, Range range, ParityError oddParity, ParityError headerParity) {
    m_packedLdf[pmt] = 0;
    if (headerParity == ERROR) m_packedLdf[pmt] |= (ERROR << HEADERPARITY_SHIFT);
    if (range == HIGH) m_packedLdf[pmt] |= (HIGH << RANGE_SHIFT); 
    if (oddParity == ERROR) m_packedLdf[pmt] |= (ERROR << ERROR_SHIFT);

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


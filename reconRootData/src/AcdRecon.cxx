#include "reconRootData/AcdRecon.h"
#include <iostream>

ClassImp(AcdRecon)

AcdRecon::AcdRecon() {
    Clear();
}

AcdRecon::~AcdRecon() {
    Clear();
}

void AcdRecon::initialize(Double_t e, Int_t count, Double_t gDoca, Double_t doca, Double_t actDist,
                          const AcdId &minDocaId, const std::vector<Double_t> &rowDoca) {
    Clear();
    m_totEnergy = e;
    m_tileCount = count;
    m_gammaDoca = gDoca;
    m_doca = doca;
    m_actDist = actDist;
    m_minDocaId = minDocaId;
    m_rowDocaCol = rowDoca;
}

void AcdRecon::initialize(Double_t e, Int_t count, Double_t gDoca, Double_t doca, Double_t actDist,
                          const AcdId &minDocaId) {
    Clear();
    m_totEnergy = e;
    m_tileCount = count;
    m_gammaDoca = gDoca;
    m_doca = doca;
    m_actDist = actDist;
    m_minDocaId = minDocaId;
}

void AcdRecon::Clear(Option_t *option) {
    m_totEnergy = 0.0;
    m_tileCount = 0;
    m_gammaDoca = -200.;
    m_doca = -200.;
    m_actDist = -200.0;
    m_rowDocaCol.clear();
}

void AcdRecon::Print(Option_t *option) const {
    TObject::Print(option);
    using namespace std;
    cout << "TileCount: " << m_tileCount;
    cout << " DOCA: " << m_doca << endl;
    m_minDocaId.Print();
    cout << " Act_Dist: " << m_actDist << endl;
    cout << "GammeDoca: " << m_gammaDoca << endl;
    cout << "RowDoca: " << endl;;
    std::vector<Double_t>::const_iterator rowDocaIt;
    for (rowDocaIt = m_rowDocaCol.begin(); rowDocaIt != m_rowDocaCol.end(); rowDocaIt++) 
        cout << (*rowDocaIt) << endl;

}
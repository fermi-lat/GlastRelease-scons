#include "reconRootData/AcdRecon.h"
#include <iostream>

ClassImp(AcdRecon)

AcdRecon::AcdRecon() {
    Clear();
}

AcdRecon::~AcdRecon() {
    Clear();
}

void AcdRecon::initialize(Double_t e, Int_t count, Double_t gDoca, 
                          Double_t doca, const AcdId &minDocaId, 
                          Double_t actDist,
                          const AcdId &maxActDistId, 
                          const std::vector<Double_t> &rowDoca,
                          const std::vector<Double_t> &rowActDist, 
			  const std::vector<AcdId> &idCol,
                          const std::vector<Double_t> &energyCol) {
    Clear();
    m_totEnergy = e;
    m_tileCount = count;
    m_gammaDoca = gDoca;
    m_doca = doca;
    m_actDist = actDist;
    m_minDocaId = minDocaId;
    m_maxActDistId = maxActDistId;
    m_rowDocaCol = rowDoca;
    m_rowActDistCol = rowActDist;
    m_idCol = idCol;
    m_energyCol = energyCol;
}

void AcdRecon::initialize(Double_t e, Int_t count, Double_t gDoca, 
                Double_t doca, const AcdId &minDocaId, Double_t actDist,
                const AcdId &maxActDistId) {
    Clear();
    m_totEnergy = e;
    m_tileCount = count;
    m_gammaDoca = gDoca;
    m_doca = doca;
    m_actDist = actDist;
    m_minDocaId = minDocaId;
    m_maxActDistId = maxActDistId;
}

void AcdRecon::Clear(Option_t* /*option*/) {
    m_totEnergy = 0.0;
    m_tileCount = 0;
    m_gammaDoca = -200.;
    m_doca = -200.;
    m_actDist = -200.0;
    m_rowDocaCol.clear();
    m_rowActDistCol.clear();
    m_idCol.clear();
    m_energyCol.clear();
}

void AcdRecon::Print(Option_t *option) const {
    TObject::Print(option);
    using namespace std;
    cout << "TileCount: " << m_tileCount;
    cout << " DOCA: " << m_doca << endl;
    m_minDocaId.Print();
    cout << " Act_Dist: " << m_actDist << endl;
    m_maxActDistId.Print();
    cout << "GammeDoca: " << m_gammaDoca << endl;
    cout << "RowDoca: " << endl;
    std::vector<Double_t>::const_iterator rowDocaIt;
    for (rowDocaIt = m_rowDocaCol.begin(); rowDocaIt != m_rowDocaCol.end(); rowDocaIt++) 
        cout << (*rowDocaIt) << endl;
    cout << "RowActDist: " << endl;
    std::vector<Double_t>::const_iterator rowActDistIt;
    for (rowActDistIt = m_rowActDistCol.begin(); rowActDistIt != m_rowActDistCol.end(); rowActDistIt++) 
        cout << (*rowActDistIt) << endl;
	cout << "Energy Collection: " << endl;
	std::vector<Double_t>::const_iterator energyIt;
	unsigned int i = 0;
	for (energyIt = m_energyCol.begin(); energyIt != m_energyCol.end(); energyIt++) {
		cout << m_idCol[i].getId() << " " << *energyIt << endl;
		i++;
	}

}

Double_t AcdRecon::getEnergy(const AcdId& id) const {
	std::vector<AcdId>::const_iterator idIt;
	unsigned int index = 0;
	for (idIt = m_idCol.begin(); idIt != m_idCol.end(); idIt++) {
		if (*idIt == id) break;
		++index;
	}
	if (idIt == m_idCol.end()) return -1.0;
	if (index >= m_energyCol.size()) return -1.0;
	return m_energyCol[index];
}

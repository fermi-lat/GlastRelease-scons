#include "reconRootData/AcdRecon.h"
#include <iostream>

ClassImp(AcdRecon)

AcdRecon::AcdRecon() {
  m_acdTkrIntersectionCol = 0;
  Clear();
}

AcdRecon::~AcdRecon() {
  Clear();
  if ( m_acdTkrIntersectionCol != 0 ) {
    m_acdTkrIntersectionCol->Delete();
  }
  delete m_acdTkrIntersectionCol;
  m_acdTkrIntersectionCol = 0;
}

void AcdRecon::initialize(Double_t e, Double_t ribbonE, Int_t count, 
                          Int_t ribbonCount, Double_t gDoca, 
                          Double_t doca, const AcdId &minDocaId, 
                          Double_t actDist,
                          const AcdId &maxActDistId, 
                          Double32_t ribbonActDist, 
                          const AcdId &ribbonActDistId,
                          const std::vector<Double_t> &rowDoca,
                          const std::vector<Double_t> &rowActDist, 
			  const std::vector<AcdId> &idCol,
                          const std::vector<Double_t> &energyCol) {
    Clear();
    m_totEnergy = e;
    m_totRibbonEnergy = ribbonE;
    m_tileCount = count;
    m_ribbonCount = ribbonCount;
    m_gammaDoca = gDoca;
    m_doca = doca;
    m_actDist = actDist;
    m_minDocaId = minDocaId;
    m_maxActDistId = maxActDistId;
    m_rowDocaCol = rowDoca;
    m_rowActDistCol = rowActDist;
    m_ribbonActDist = ribbonActDist;
    m_ribbonActDistId = ribbonActDistId;
    m_idCol = idCol;
    m_energyCol = energyCol;
}

void AcdRecon::initialize(Double_t e, Double_t ribbonE, Int_t count, 
                Int_t ribbonCount, Double_t gDoca, 
                Double_t doca, const AcdId &minDocaId, Double_t actDist,
                const AcdId &maxActDistId) {
    Clear();
    m_totEnergy = e;
    m_totRibbonEnergy = ribbonE;
    m_tileCount = count;
    m_ribbonCount = ribbonCount;
    m_gammaDoca = gDoca;
    m_doca = doca;
    m_actDist = actDist;
    m_minDocaId = minDocaId;
    m_maxActDistId = maxActDistId;
}

void AcdRecon::Clear(Option_t* /*option*/) {
    m_totEnergy = 0.0;
    m_totRibbonEnergy = 0.0;
    m_tileCount = 0;
    m_ribbonCount = 0;
    m_gammaDoca = -200.;
    m_doca = -200.;
    m_actDist = -200.0;
    m_rowDocaCol.clear();
    m_rowActDistCol.clear();
    m_idCol.clear();
    m_energyCol.clear();
    if ( m_acdTkrIntersectionCol ) {
      m_acdTkrIntersectionCol->Clear();
    }
}

void AcdRecon::Print(Option_t *option) const {
    TObject::Print(option);
    using namespace std;
    cout << "TileCount: " << m_tileCount;
    cout << "RibbonCount: " << m_ribbonCount;
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


Bool_t AcdRecon::hasHit(UInt_t acdId) const {  
  std::vector<AcdId>::const_iterator idIt;
  for (idIt = m_idCol.begin(); idIt != m_idCol.end(); idIt++) {
    if ( idIt->getId() == acdId ) return kTRUE;
  }
  return kFALSE;
}


AcdTkrIntersection* AcdRecon::addAcdTkrIntersection(AcdTkrIntersection& inter) {
  UInt_t numInter = nAcdIntersections();
  if ( m_acdTkrIntersectionCol == 0 ) {
    m_acdTkrIntersectionCol = new TClonesArray( AcdTkrIntersection::Class_Name() );
  }
  // fill the TClonesArray w/ a placement new call
  TClonesArray& intersects = *m_acdTkrIntersectionCol;
  new (intersects[numInter]) AcdTkrIntersection(inter);  
  return (AcdTkrIntersection*)(intersects[numInter]);
}


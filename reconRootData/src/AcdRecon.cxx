#include "reconRootData/AcdRecon.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(AcdRecon)

static const Double_t maxDoca = 2000.0;

AcdRecon::AcdRecon() {
  m_acdTkrIntersectionCol = 0;
  m_acdTkrPocaCol = 0;
  m_acdHitCol = 0;
  Clear();
}

AcdRecon::~AcdRecon() {
  Clear();
  if ( m_acdTkrIntersectionCol != 0 ) {
    m_acdTkrIntersectionCol->Delete();
  }
  delete m_acdTkrIntersectionCol;
  m_acdTkrIntersectionCol = 0;
  if ( m_acdTkrPocaCol != 0 ) {
    m_acdTkrPocaCol->Delete();
  }
  delete m_acdTkrPocaCol;
  m_acdTkrPocaCol = 0;
  if ( m_acdHitCol != 0 ) {
    m_acdHitCol->Delete();
  }
  delete m_acdHitCol;
  m_acdHitCol = 0;

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
                          const std::vector<Double_t> &energyCol,
                          Double32_t cornerDoca) {
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
    m_cornerDoca = cornerDoca;
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
    m_gammaDoca = maxDoca;
    m_doca = maxDoca;
    m_actDist = -maxDoca;
    m_rowDocaCol.clear();
    m_rowActDistCol.clear();
    m_idCol.clear();
    m_energyCol.clear();
    if ( m_acdTkrIntersectionCol ) {
      m_acdTkrIntersectionCol->Clear();
    }
    if ( m_acdTkrPocaCol ) {
      m_acdTkrPocaCol->Clear();
    }
    if ( m_acdHitCol ) {
      m_acdHitCol->Clear();
    }
     m_ribbonActDist = -maxDoca;
    m_cornerDoca = maxDoca;
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


AcdTkrPoca* AcdRecon::addAcdTkrPoca(AcdTkrPoca& poca) {
  UInt_t numPoca = nAcdTkrPoca();
  if ( m_acdTkrPocaCol == 0 ) {
    m_acdTkrPocaCol = new TClonesArray( AcdTkrPoca::Class_Name() );
  }
  // fill the TClonesArray w/ a placement new call
  TClonesArray& pocas = *m_acdTkrPocaCol;
  new (pocas[numPoca]) AcdTkrPoca(poca);
  return (AcdTkrPoca*)(pocas[numPoca]);
}
    

AcdHit* AcdRecon::addAcdHit(AcdHit& hit) {
  UInt_t numHit = nAcdHit();
  if ( m_acdHitCol == 0 ) {
    m_acdHitCol = new TClonesArray( AcdHit::Class_Name() );
  }
  // fill the TClonesArray w/ a placement new call
  TClonesArray& hits = *m_acdHitCol;
  new (hits[numHit]) AcdHit(hit);
  return (AcdHit*)(hits[numHit]);
}
    

//======================================================
// For Unit Tests
//======================================================

void AcdRecon::Fake( Int_t ievent, Float_t randNum ) {

    Float_t f = Float_t(ievent);
    Double_t energy = f ;
    Double_t ribbonE = energy;
    Int_t count = 5;
    Int_t ribbonCount = 2;
    Double_t gDoca = f*randNum;
    Double_t doca = randNum;
    Double_t actDist = 2.*f;
    Double_t ribActDist = 3.*f;
    AcdId minDocaId(0, 0, 3, 2);
    AcdId maxActDistId(0,0,2,2);
    AcdId ribActDistId(0,1,4,3);
    ribActDistId.setReadoutIndex(ievent) ;
    std::vector<Double_t> rowDocaCol;
    rowDocaCol.push_back(randNum);
    rowDocaCol.push_back(f);
    std::vector<Double_t> rowActDistCol ;
    rowActDistCol.push_back(randNum) ;
    rowActDistCol.push_back(f) ;
    std::vector<AcdId> idCol ;
    idCol.push_back(AcdId(0, 0, 3, 2)) ;
    idCol.push_back(AcdId(0, 2, 2, 1)) ;
    std::vector<Double_t> energyCol ;
    energyCol.push_back(f) ;
    energyCol.push_back(2.*f) ;
    Double32_t cornerDoca = f*randNum ;
    
    initialize(
        energy, ribbonE, count, ribbonCount, gDoca, doca, 
        minDocaId, actDist, maxActDistId,
        ribActDist, ribActDistId,
        rowDocaCol, rowActDistCol,
        idCol, energyCol, cornerDoca
    ) ;
 }

#define UTIL_COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),a.get ## att(),#att)
#define ID_COMPARE_IN_RANGE(att) get ## att().CompareInRange(a.get ## att(),#att)

Bool_t AcdRecon::CompareInRange( const AcdRecon & a, const std::string & name ) const {

    bool result = true ;

    result = UTIL_COMPARE_IN_RANGE(Energy) && result ;
    result = UTIL_COMPARE_IN_RANGE(RibbonEnergy) && result ;
    result = UTIL_COMPARE_IN_RANGE(TileCount) && result ;
    result = UTIL_COMPARE_IN_RANGE(RibbonCount) && result ;
    result = UTIL_COMPARE_IN_RANGE(GammaDoca) && result ;
    result = UTIL_COMPARE_IN_RANGE(Doca) && result ;
    result = ID_COMPARE_IN_RANGE(MinDocaId) && result ;
    result = UTIL_COMPARE_IN_RANGE(ActiveDist) && result ;
    result = UTIL_COMPARE_IN_RANGE(RibbonActiveDist) && result ;
    result = ID_COMPARE_IN_RANGE(MinDocaId) && result ;
    result = ID_COMPARE_IN_RANGE(MaxActDistId) && result ;
    result = ID_COMPARE_IN_RANGE(RibbonActDistId) && result ;
    result = UTIL_COMPARE_IN_RANGE(RowDocaCol) && result ;
    result = UTIL_COMPARE_IN_RANGE(RowActDistCol) && result ;
    result = UTIL_COMPARE_IN_RANGE(EnergyCol) && result ;
    result = UTIL_COMPARE_IN_RANGE(IdCol) && result ;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }

    return result ;

}


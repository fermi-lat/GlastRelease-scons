
#include "reconRootData/AcdTkrAssoc.h"
#include "reconRootData/AcdTkrHitPocaV2.h"
#include "reconRootData/AcdTkrGapPocaV2.h"

#include "Riostream.h"


ClassImp(AcdTkrAssoc) ;

/// Default constructor.  Set everything to default values
AcdTkrAssoc::AcdTkrAssoc() 
  :TObject(),
   m_hitPocae(AcdTkrHitPocaV2::Class()),
   m_gapPocae(AcdTkrGapPocaV2::Class()){   
  Clear("");
}

/// Copy constructor
AcdTkrAssoc::AcdTkrAssoc(const AcdTkrAssoc& other)
  :TObject(other),
   m_trackIndex(other.m_trackIndex),
   m_upward(other.m_upward),
   m_energy(other.m_energy),
   m_start(other.m_start),
   m_dir(other.m_dir),
   m_arcLength(other.m_arcLength),
   m_cov_start(other.m_cov_start),
   m_cov_end(other.m_cov_end),
   m_tkrSSDVeto(other.m_tkrSSDVeto),
   m_cornerDoca(other.m_cornerDoca),
   m_hitPocae(AcdTkrHitPocaV2::Class()),
   m_gapPocae(AcdTkrGapPocaV2::Class()){
}

/// Assignment operator
AcdTkrAssoc& AcdTkrAssoc::operator=(const AcdTkrAssoc& other)
{
  if ( this == &other ) return *this;
  m_trackIndex = other.m_trackIndex;
  m_upward = other.m_upward;
  m_energy = other.m_energy;
  m_start = other.m_start;
  m_dir = other.m_dir;
  m_arcLength = other.m_arcLength;
  m_cov_start = other.m_cov_start;
  m_cov_end = other.m_cov_end;
  m_tkrSSDVeto = other.m_tkrSSDVeto;
  m_cornerDoca = other.m_cornerDoca;
  return *this;
}

const AcdTkrHitPocaV2* AcdTkrAssoc::getHitPoca(UInt_t i) const {
  if ( i >= (UInt_t)m_hitPocae.GetEntriesFast() ) return 0;
  const TObject* obj = m_hitPocae.UncheckedAt(i);
  const AcdTkrHitPocaV2* ret = (const AcdTkrHitPocaV2*)(obj);
  return ret;
}

const AcdTkrGapPocaV2* AcdTkrAssoc::getGapPoca(UInt_t i) const {
  if ( i >= (UInt_t)m_gapPocae.GetEntriesFast() ) return 0;
  const TObject* obj = m_gapPocae.UncheckedAt(i);
  const AcdTkrGapPocaV2* ret = (const AcdTkrGapPocaV2*)(obj);
  return ret;
}


AcdTkrHitPocaV2* AcdTkrAssoc::addHitPoca(const AcdTkrHitPocaV2& poca) {
  UInt_t n = m_hitPocae.GetEntriesFast();
  AcdTkrHitPocaV2* newPoca = (AcdTkrHitPocaV2*)(new (m_hitPocae[n]) AcdTkrHitPocaV2(poca));
  return newPoca;
}

AcdTkrGapPocaV2* AcdTkrAssoc::addGapPoca(const AcdTkrGapPocaV2& poca) {
  UInt_t n = m_gapPocae.GetEntriesFast();
  AcdTkrGapPocaV2* newPoca = (AcdTkrGapPocaV2*)(new (m_gapPocae[n]) AcdTkrGapPocaV2(poca));
  return newPoca;
}

void AcdTkrAssoc::setPoint(const AcdTkrPointV2& point) {
  m_point = point;
}

/// set all the values
void AcdTkrAssoc::set(Int_t index, Bool_t up, Float_t energy, 
                      const TVector3& start, const TVector3& dir, Float_t arcLength,
                      const TMatrixDSym& covStart, const TMatrixDSym& covEnd,
                      Int_t tkrSSDVeto, Float_t cornerDoca)
{
  m_trackIndex = index;
  m_upward = up;
  m_energy = energy;
  m_start = start;
  m_dir = dir;
  m_arcLength = arcLength;
  m_cov_start = covStart;
  m_cov_end = covEnd;
  m_tkrSSDVeto = tkrSSDVeto;
  m_cornerDoca = cornerDoca;
}  


/// reset all the values to their default
void AcdTkrAssoc::Clear(Option_t* option) 
{
  TObject::Clear(option);
  m_trackIndex = -1;  
  m_energy = 0.;;
  m_start.Clear(option);
  m_dir.Clear(option);
  m_arcLength = 0.;
  m_cov_start.Clear(option);
  m_cov_end.Clear(option);
  m_tkrSSDVeto = 0;
  m_cornerDoca = 0.;  
}

/// Print out this structure
void AcdTkrAssoc::Print(Option_t *option) const 
{
  std::cout << "AcdTkrAssoc.  " << m_trackIndex << ' ' << (m_upward ? "up" : "down") << " E=" << m_energy 
            << " s= " << m_arcLength << " SSDVeto = " << m_tkrSSDVeto 
            << " CornerDoca = " << m_cornerDoca << std::endl;
  m_start.Print(option);
  m_dir.Print(option);
  m_cov_start.Print(option);
  m_cov_end.Print(option);
}

void AcdTkrAssoc::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */) {
}

Bool_t AcdTkrAssoc::CompareInRange( const AcdTkrAssoc&, const std::string & /* name */ ) const {
  return kTRUE;
}

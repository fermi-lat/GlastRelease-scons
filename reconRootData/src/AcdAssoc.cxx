
#include "reconRootData/AcdAssoc.h"
#include "reconRootData/AcdTkrHitPocaV2.h"
#include "reconRootData/AcdTkrGapPocaV2.h"

#include "Riostream.h"


ClassImp(AcdAssoc) ;

/// Default constructor.  Set everything to default values
AcdAssoc::AcdAssoc() 
  :TObject(),
   m_trackIndex(-1),
   m_upward(kFALSE),
   m_energy(0.),
   m_arcLength(0.),
   m_cov_start(5),
   m_cov_end(5),
   m_tkrSSDVeto(-1),
   m_cornerDoca(0.),   
   m_hitPocae(AcdTkrHitPocaV2::Class()),
   m_gapPocae(AcdTkrGapPocaV2::Class()),
   m_energy15(0.),
   m_energy30(0.),
   m_energy45(0.),
   m_triggerEnergy15(0.),
   m_triggerEnergy30(0.),
   m_triggerEnergy45(0.){   
  Clear("");
}

/// Copy constructor
AcdAssoc::AcdAssoc(const AcdAssoc& other)
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
   m_gapPocae(AcdTkrGapPocaV2::Class()),
   m_point(other.m_point),
   m_energy15(other.m_energy15),
   m_energy30(other.m_energy30),
   m_energy45(other.m_energy45),
   m_triggerEnergy15(other.m_triggerEnergy15),
   m_triggerEnergy30(other.m_triggerEnergy30),
   m_triggerEnergy45(other.m_triggerEnergy45){      
}


const AcdTkrHitPocaV2* AcdAssoc::getHitPoca(UInt_t i) const {
  if ( i >= (UInt_t)m_hitPocae.GetEntriesFast() ) return 0;
  const TObject* obj = m_hitPocae.UncheckedAt(i);
  const AcdTkrHitPocaV2* ret = (const AcdTkrHitPocaV2*)(obj);
  return ret;
}

const AcdTkrGapPocaV2* AcdAssoc::getGapPoca(UInt_t i) const {
  if ( i >= (UInt_t)m_gapPocae.GetEntriesFast() ) return 0;
  const TObject* obj = m_gapPocae.UncheckedAt(i);
  const AcdTkrGapPocaV2* ret = (const AcdTkrGapPocaV2*)(obj);
  return ret;
}


AcdTkrHitPocaV2* AcdAssoc::addHitPoca(const AcdTkrHitPocaV2& poca) {
  UInt_t n = m_hitPocae.GetEntriesFast();
  AcdTkrHitPocaV2* newPoca = (AcdTkrHitPocaV2*)(new (m_hitPocae[n]) AcdTkrHitPocaV2(poca));
  return newPoca;
}

AcdTkrGapPocaV2* AcdAssoc::addGapPoca(const AcdTkrGapPocaV2& poca) {
  UInt_t n = m_gapPocae.GetEntriesFast();
  AcdTkrGapPocaV2* newPoca = (AcdTkrGapPocaV2*)(new (m_gapPocae[n]) AcdTkrGapPocaV2(poca));
  return newPoca;
}

void AcdAssoc::setPoint(const AcdTkrPointV2& point) {
  m_point = point;
}

/// set all the values
void AcdAssoc::set(Int_t index, Bool_t up, Float_t energy, 
		   const TVector3& start, const TVector3& dir, Float_t arcLength,
		   const TMatrixFSym& covStart, const TMatrixFSym& covEnd,
 		   Int_t tkrSSDVeto, Float_t cornerDoca,
		   float energy15, float energy30, float energy45,
		   float triggerEnergy15, float triggerEnergy30, float triggerEnergy45)
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
  m_energy15 = energy15;
  m_energy30 = energy30;
  m_energy45 = energy45;
  m_triggerEnergy15 = triggerEnergy15;
  m_triggerEnergy30 = triggerEnergy30;
  m_triggerEnergy45 = triggerEnergy45;      
}  


/// reset all the values to their default
void AcdAssoc::Clear(Option_t* option) 
{
  TObject::Clear(option);
  m_trackIndex = -1;  
  m_energy = 0.;;
  m_start.Clear(option);
  m_dir.Clear(option);
  m_arcLength = 0.;
  m_cov_start.Zero();
  m_cov_end.Zero();
  m_tkrSSDVeto = 0;
  m_cornerDoca = 0.;  
  m_energy15 = 0.;
  m_energy30 = 0.;
  m_energy45 = 0.;
  m_triggerEnergy15 = 0.;
  m_triggerEnergy30 = 0.;
  m_triggerEnergy45 = 0.;      
}

/// Print out this structure
void AcdAssoc::Print(Option_t *option) const 
{
  std::cout << "AcdAssoc.  " << m_trackIndex << ' ' << (m_upward ? "up" : "down") << " E=" << m_energy 
            << " s= " << m_arcLength << " SSDVeto = " << m_tkrSSDVeto 
            << " CornerDoca = " << m_cornerDoca << std::endl;
  m_start.Print(option);
  m_dir.Print(option);
  m_cov_start.Print(option);
  m_cov_end.Print(option);
}

void AcdAssoc::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */) {
}

Bool_t AcdAssoc::CompareInRange( const AcdAssoc&, const std::string & /* name */ ) const {
  return kTRUE;
}

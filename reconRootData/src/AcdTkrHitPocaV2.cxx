#include "reconRootData/AcdTkrHitPocaV2.h"
#include "Riostream.h"

ClassImp(AcdTkrHitPocaV2) ;

/// Default constructor.  Set everything to default values
AcdTkrHitPocaV2::AcdTkrHitPocaV2() 
  :TObject(),AcdTkrLocalCoordsV2(),AcdPocaDataV2(){
  Clear("");
}

/// Copy constructor
AcdTkrHitPocaV2::AcdTkrHitPocaV2(const AcdTkrHitPocaV2& other)
  :TObject(other),AcdTkrLocalCoordsV2(other),AcdPocaDataV2(other)
{
  //Clear();
  set(other.getId(),other.getTrackIndex(),other.m_mips,
      other.vetoSigmaHit(),other.vetoSigmaProj(),other.vetoSigmaProp());
}

/// Assignment operator
AcdTkrHitPocaV2& AcdTkrHitPocaV2::operator=(const AcdTkrHitPocaV2& other)
{
  if ( this == &other ) return *this;
  set(other.getId(),other.getTrackIndex(),other.m_mips,
      other.vetoSigmaHit(),other.vetoSigmaProj(),other.vetoSigmaProp());
  AcdTkrLocalCoordsV2::copy(other);
  AcdPocaDataV2::copy(other);
  return *this;
}


/// set all the values
void AcdTkrHitPocaV2::set(const AcdId& acdId, int trackIndex,
                        const float mips[2],
                        Float_t vetoSigmaHit, Float_t vetoSigmaProj, Float_t vetoSigmaProp)
{
  m_id = acdId;
  m_trackIndex = trackIndex;
  m_mips[0] = mips[0];
  m_mips[1] = mips[1];
  m_vetoSigmaHit = vetoSigmaHit;
  m_vetoSigmaProj = vetoSigmaProj;
  m_vetoSigmaProp = vetoSigmaProp;
}  


/// reset all the values to their default
void AcdTkrHitPocaV2::Clear(Option_t* option) 
{
  TObject::Clear(option);
  m_id = AcdId();
  m_trackIndex = -1;  
  m_mips[0] = 0.;
  m_mips[1] = 0.;
  m_vetoSigmaHit = 0.;
  m_vetoSigmaProj = 0.;
  m_vetoSigmaProp = 0.;
  AcdTkrLocalCoordsV2::clear();
  AcdPocaDataV2::clear();
}

/// Print out this structure
void AcdTkrHitPocaV2::Print(Option_t *option) const 
{
  using namespace std;    
  cout << "AcdTkrHitPocaV2.  id: " << m_id.getId() << "  tk: " << m_trackIndex  
       << " Mips("  << m_mips[0] << ',' << m_mips[1] << ')' << endl;
  AcdTkrLocalCoordsV2::print(option);
  AcdPocaDataV2::print(option);
}

void AcdTkrHitPocaV2::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */) {
}

Bool_t AcdTkrHitPocaV2::CompareInRange( const AcdTkrHitPocaV2&, const std::string & /* name */ ) const {
  return kTRUE;
}

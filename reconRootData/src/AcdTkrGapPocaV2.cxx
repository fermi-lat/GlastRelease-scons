#include "reconRootData/AcdTkrGapPocaV2.h"
#include "Riostream.h"

ClassImp(AcdTkrGapPocaV2) ;

/// Default constructor.  Set everything to default values
AcdTkrGapPocaV2::AcdTkrGapPocaV2() 
  :TObject(),AcdTkrLocalCoordsV2(),AcdPocaDataV2(){
  Clear("");
}


/// Copy constructor
AcdTkrGapPocaV2::AcdTkrGapPocaV2(const AcdTkrGapPocaV2& other)
  :TObject(other),AcdTkrLocalCoordsV2(other),AcdPocaDataV2(other)
{
  //Clear();
  initialize(other.getGapId(),other.getTrackIndex(),
             other.vetoSigmaHit(),other.vetoSigmaProj(),other.vetoSigmaProp());
}

/// Assignment operator
AcdTkrGapPocaV2& AcdTkrGapPocaV2::operator=(const AcdTkrGapPocaV2& other)
{
  if ( this == &other ) return *this;
  initialize(other.getGapId(),other.getTrackIndex(),
             other.vetoSigmaHit(),other.vetoSigmaProj(),other.vetoSigmaProp());
  AcdTkrLocalCoordsV2::copy(other);
  AcdPocaDataV2::copy(other); 
  return *this;
}


/// set all the values
void AcdTkrGapPocaV2::initialize(const AcdGapId& acdId, Int_t trackIndex,
                               Float_t vetoSigmaHit, Float_t vetoSigmaProj, Float_t vetoSigmaProp)
{
  m_id = acdId;
  m_trackIndex = trackIndex;
  m_vetoSigmaHit = vetoSigmaHit;
  m_vetoSigmaProj = vetoSigmaProj;
  m_vetoSigmaProp = vetoSigmaProp;
}  


/// reset all the values to their default
void AcdTkrGapPocaV2::Clear(Option_t* option) 
{
  TObject::Clear(option);
  m_id = AcdGapId();
  m_trackIndex = -1;
  m_vetoSigmaHit = 0.;
  m_vetoSigmaProj = 0.;
  m_vetoSigmaProp = 0.;
  AcdTkrLocalCoordsV2::clear();
  AcdPocaDataV2::clear();
}

/// Print out this structure
void AcdTkrGapPocaV2::Print(Option_t *option) const 
{
  using namespace std;    
  cout << "AcdTkrGapPocaV2.  id: " << m_id.asDecimal() << "  tk: " << m_trackIndex << endl;
  AcdTkrLocalCoordsV2::print(option);
  AcdPocaDataV2::print(option);
}

void AcdTkrGapPocaV2::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */) {
}

Bool_t AcdTkrGapPocaV2::CompareInRange( const AcdTkrGapPocaV2&, const std::string & /* name */ ) const {
  return kTRUE;
}

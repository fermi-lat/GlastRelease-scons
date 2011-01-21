#include "reconRootData/AcdTkrHitPoca.h"
#include "Riostream.h"

ClassImp(AcdTkrHitPoca) ;

/// Default constructor.  Set everything to default values
AcdTkrHitPoca::AcdTkrHitPoca() 
  :TObject(),AcdTkrLocalCoords(),AcdPocaData(){
  Clear("");
}

/// Constructor for use in transient -> persistent conversion 
/// Takes arguements as they are stored in ROOT
AcdTkrHitPoca::AcdTkrHitPoca(const AcdId& acdId, int trackIndex,
                             const AcdTkrLocalCoords& local,
                             const AcdPocaData& pocaData)
  :TObject(),AcdTkrLocalCoords(),AcdPocaData()
{
  //Clear("");
  set(acdId,trackIndex,local,pocaData);
}

/// Copy constructor
AcdTkrHitPoca::AcdTkrHitPoca(const AcdTkrHitPoca& other)
  :TObject(other),AcdTkrLocalCoords(),AcdPocaData()
{
  //Clear();
  set(other.getId(),other.getTrackIndex(),other,other);
}

/// Assignment operator
AcdTkrHitPoca& AcdTkrHitPoca::operator=(const AcdTkrHitPoca& other)
{
  if ( this == &other ) return *this;
  set(other.getId(),other.getTrackIndex(),other,other);
  return *this;
}


/// set all the values
void AcdTkrHitPoca::set(const AcdId& acdId, Int_t trackIndex,
                        const AcdTkrLocalCoords& local,
                        const AcdPocaData& pocaData)
{
  m_id = acdId;
  m_trackIndex = trackIndex;
  AcdTkrLocalCoords::copy(local);
  AcdPocaData::copy(pocaData);
}  


/// reset all the values to their default
void AcdTkrHitPoca::Clear(Option_t* option) 
{
  TObject::Clear(option);
  m_id = AcdId();
  m_trackIndex = -1;
  AcdTkrLocalCoords::clear();
  AcdPocaData::clear();
}

/// Print out this structure
void AcdTkrHitPoca::Print(Option_t *option) const 
{
  using namespace std;    
  cout << "AcdTkrHitPoca.  id: " << m_id.getId() << "  tk: " << m_trackIndex << endl;;
  AcdTkrLocalCoords::print(option);
  AcdPocaData::print(option);
}

void AcdTkrHitPoca::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */) {
}

Bool_t AcdTkrHitPoca::CompareInRange( const AcdTkrHitPoca&, const std::string & /*name*/ ) const {
  return kTRUE;
}

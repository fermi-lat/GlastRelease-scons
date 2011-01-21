#include "reconRootData/AcdTkrGapPoca.h"
#include "Riostream.h"

ClassImp(AcdTkrGapPoca) ;

/// Default constructor.  Set everything to default values
AcdTkrGapPoca::AcdTkrGapPoca() 
  :TObject(),AcdTkrLocalCoords(),AcdPocaData(){
  Clear("");
}

/// Constructor for use in transient -> persistent conversion 
/// Takes arguements as they are stored in ROOT
AcdTkrGapPoca::AcdTkrGapPoca(const AcdGapId& acdId, int trackIndex,
                             const AcdTkrLocalCoords& local,
                             const AcdPocaData& pocaData)
  :TObject(),AcdTkrLocalCoords(),AcdPocaData()
{
  //Clear("");
  initialize(acdId,trackIndex,local,pocaData);
}

/// Copy constructor
AcdTkrGapPoca::AcdTkrGapPoca(const AcdTkrGapPoca& other)
  :TObject(other),AcdTkrLocalCoords(),AcdPocaData()
{
  //Clear();
  initialize(other.getGapId(),other.getTrackIndex(),other,other);
}

/// Assignment operator
AcdTkrGapPoca& AcdTkrGapPoca::operator=(const AcdTkrGapPoca& other)
{
  if ( this == &other ) return *this;
  initialize(other.getGapId(),other.getTrackIndex(),other,other);
  return *this;
}


/// set all the values
void AcdTkrGapPoca::initialize(const AcdGapId& acdId, Int_t trackIndex,
                               const AcdTkrLocalCoords& local,
                               const AcdPocaData& pocaData)
{
  m_id = acdId;
  m_trackIndex = trackIndex;
  AcdTkrLocalCoords::copy(local);
  AcdPocaData::copy(pocaData);
}  


/// reset all the values to their default
void AcdTkrGapPoca::Clear(Option_t* option) 
{
  TObject::Clear(option);
  m_id = AcdGapId();
  m_trackIndex = -1;
  AcdTkrLocalCoords::clear();
  AcdPocaData::clear();
}

/// Print out this structure
void AcdTkrGapPoca::Print(Option_t *option) const 
{
  using namespace std;    
  cout << "AcdTkrGapPoca.  id: " << m_id.asDecimal() << "  tk: " << m_trackIndex << endl;
  AcdTkrLocalCoords::print(option);
  AcdPocaData::print(option);
}

void AcdTkrGapPoca::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */) {
}

Bool_t AcdTkrGapPoca::CompareInRange( const AcdTkrGapPoca&, const std::string & /*name*/ ) const {
  return kTRUE;
}

#include "reconRootData/AcdTkrPointV2.h"
#include "Riostream.h"

ClassImp(AcdTkrPointV2) ;

/// Default constructor.  Set everything to default values
AcdTkrPointV2::AcdTkrPointV2() {
  Clear("");
}

/// Copy constructor
AcdTkrPointV2::AcdTkrPointV2(const AcdTkrPointV2& other)
  :TObject(other),
   AcdTkrLocalCoordsV2(other){
  //Clear();
  set(other.getTrackIndex());
}

/// Assignment operator
AcdTkrPointV2& AcdTkrPointV2::operator=(const AcdTkrPointV2& other)
{
  if ( this == &other ) return *this;
  set(other.getTrackIndex());
  AcdTkrLocalCoordsV2::copy(other);
  return *this;
}


/// set all the values
void AcdTkrPointV2::set(Int_t trackIndex)			    
{
  m_trackIndex = trackIndex;
}  

void AcdTkrPointV2::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */) {
}

Bool_t AcdTkrPointV2::CompareInRange( const AcdTkrPointV2&, const std::string & /* name */) const 
{
  // FIXME
  return kTRUE;
}

/// reset all the values to their default
void AcdTkrPointV2::Clear(Option_t* /* option */) 
{
  m_trackIndex = -1;
  AcdTkrLocalCoordsV2::clear();
}

/// Print out this structure
void AcdTkrPointV2::Print(Option_t *option) const 
{
  std::cout << "AcdTkrPointV2.  " << m_trackIndex << std::endl;
  AcdTkrLocalCoordsV2::print(option);
}


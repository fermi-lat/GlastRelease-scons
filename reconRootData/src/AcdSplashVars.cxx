#include "reconRootData/AcdSplashVars.h"
#include "Riostream.h"

ClassImp(AcdSplashVars) ;

/// Default constructor.  Set everything to default values
AcdSplashVars::AcdSplashVars() 
  :TObject(){
  Clear("");
}

/// Constructor for use in transient -> persistent conversion 
/// Takes arguements as they are stored in ROOT
AcdSplashVars::AcdSplashVars(const AcdId& acdId, int trackIndex, 
			 const TVector3& calEntryPoint, const TVector3& calEntryVector,
			 const Float_t& tileSolidAngle, const Float_t& weightedTrackAngle,
			 const Float_t& weightedPathlength)
  :TObject()
{
  //Clear("");
  set(acdId,trackIndex,
      calEntryPoint,calEntryVector,
      tileSolidAngle,weightedTrackAngle,weightedPathlength);
}

/// Copy constructor
AcdSplashVars::AcdSplashVars(const AcdSplashVars& other)
  :TObject(other)
{
  //Clear();
  set(other.getId(),other.getTrackIndex(),
      other.calEntryPoint(),other.calEntryVector(),
      other.tileSolidAngle(),other.weightedTrackAngle(),other.weightedPathlength());
}

/// Assignment operator
AcdSplashVars& AcdSplashVars::operator=(const AcdSplashVars& other)
{
  if ( this == &other ) return *this;
  set(other.getId(),other.getTrackIndex(),
      other.calEntryPoint(),other.calEntryVector(),
      other.tileSolidAngle(),other.weightedTrackAngle(),other.weightedPathlength());
  return *this;
}


/// set all the values
void AcdSplashVars::set(const AcdId& acdId, Int_t trackIndex,
		      const TVector3& calEntryPoint, const TVector3& calEntryVector,
		      const Float_t& tileSolidAngle, const Float_t& weightedTrackAngle,
		      const Float_t& weightedPathlength)
{
  m_id = acdId;
  m_trackIndex = trackIndex;

  m_calEntryPoint = calEntryPoint;
  m_calEntryVector = calEntryVector;
  m_tileSolidAngle = tileSolidAngle;
  m_weightedTrackAngle = weightedTrackAngle;
  m_weightedPathlength = weightedPathlength;
}  


/// reset all the values to their default
void AcdSplashVars::Clear(Option_t* option) 
{
  TObject::Clear(option);
  m_id = AcdId();
  m_trackIndex = -1;

  m_calEntryPoint.Clear();
  m_calEntryVector.Clear();
  m_tileSolidAngle = -1.;
  m_weightedTrackAngle = -99.;
  m_weightedPathlength = -1.;
}

/// Print out this structure
void AcdSplashVars::Print(Option_t *option) const 
{
  using namespace std;    
  cout << "AcdSplashVars.  id: " << m_id.getId() << "  tk: " << m_trackIndex << endl;;
  
  m_calEntryPoint.Print(option);
  m_calEntryVector.Print(option);

  cout << "Solid Angle: " << m_tileSolidAngle << endl;
  cout << "Weigted Track Angle: " << m_weightedTrackAngle << endl;
  cout << "Weigted Path Length: " << m_weightedPathlength << endl;
}

void AcdSplashVars::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */) {
}

Bool_t AcdSplashVars::CompareInRange( const AcdSplashVars&, const std::string & name ) const {
  return kTRUE;
}

#include "reconRootData/AcdTkrPoint.h"
#include "Riostream.h"

ClassImp(AcdTkrPoint) ;

/// Default constructor.  Set everything to default values
AcdTkrPoint::AcdTkrPoint() {
  Clear("");
}

/// Constructor for use in transient -> persistent conversion 
/// Takes arguements as they are stored in ROOT
AcdTkrPoint::AcdTkrPoint(double arcLength, int face, 
                         const TVector3& point, const TkrTrackParams& paramsAtPoint)
{
  //Clear();
  set(arcLength,face,point,paramsAtPoint);
}

/// Copy constructor
AcdTkrPoint::AcdTkrPoint(const AcdTkrPoint& other)
  :TObject(other)
{
  //Clear();
  set(other.getArcLength(),other.getFace(),
      other.getPoint(),other.getParamsAtPoint());
}

/// Assignment operator
AcdTkrPoint& AcdTkrPoint::operator=(const AcdTkrPoint& other)
{
  if ( this == &other ) return *this;
  set(other.getArcLength(),other.getFace(),
      other.getPoint(),other.getParamsAtPoint());
  return *this;
}


/// set all the values
void AcdTkrPoint::set(double arcLength, int face, 
                      const TVector3& point, const TkrTrackParams& paramsAtPoint)                            
{
  m_arcLength = arcLength;
  m_face = face;
  m_point = point;
  m_paramsAtPoint = paramsAtPoint;
}  

void AcdTkrPoint::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */) {
}

Bool_t AcdTkrPoint::CompareInRange( const AcdTkrPoint&, const std::string & /* name */) const 
{
  // FIXME
  return kTRUE;
}

/// reset all the values to their default
void AcdTkrPoint::Clear(Option_t* option) 
{
  m_arcLength = 0.;
  m_face = 0;
  m_point.SetXYZ(0.,0.,0.);
  m_paramsAtPoint.Clear(option);
}

/// Print out this structure
void AcdTkrPoint::Print(Option_t *option) const 
{
  using namespace std;    
  cout << "AcdTkrPoint.  s: " << m_arcLength << "  face: " << (int)m_face
       << "  x: [" << m_point.X() << ' ' << m_point.Y() << ' ' << m_point.Z() << "] " << endl;
  m_paramsAtPoint.Print(option);
}


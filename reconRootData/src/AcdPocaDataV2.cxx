#include "reconRootData/AcdPocaDataV2.h"
#include "Riostream.h"

ClassImp(AcdPocaDataV2) ;

/// Default constructor.  Set everything to default values
AcdPocaDataV2::AcdPocaDataV2() {
  clear();
}

/// Constructor for use in transient -> persistent conversion 
/// Takes arguements as they are stored in ROOT
AcdPocaDataV2::AcdPocaDataV2(Int_t volume, Int_t region, Float_t arcLength, 
                             Float_t doca, Float_t docaErrProj,  Float_t docaErrProp,
                             const TVector3& poca, const TVector3& pocaVector){
  //Clear();
  set(volume,region,arcLength,doca,docaErrProj,docaErrProp,poca,pocaVector);
}

/// Copy constructor
AcdPocaDataV2::AcdPocaDataV2(const AcdPocaDataV2& other)
{
  //Clear();
  set(other.getVolume(),other.getRegion(),other.getArcLength(),
      other.getDoca(),other.getDocaErrProj(),other.getDocaErrProp(),
      other.getPoca(),other.getPocaVector());             
}

/// Assignment operator
AcdPocaDataV2& AcdPocaDataV2::operator=(const AcdPocaDataV2& other)
{
  if ( this == &other ) return *this;
  set(other.getVolume(),other.getRegion(),other.getArcLength(),
      other.getDoca(),other.getDocaErrProj(),other.getDocaErrProp(),
      other.getPoca(),other.getPocaVector());         
  return *this;
}


/// set all the values
void AcdPocaDataV2::set(Int_t volume, Int_t region, Float_t arcLength,
                      Float_t doca, Float_t docaErrProj, Float_t docaErrProp,
                      const TVector3& poca, const TVector3& pocaVector)                       
{
  m_volume = volume;
  m_region = region;
  m_arcLength = arcLength;
  m_doca = doca;
  m_docaErr_proj = docaErrProj;
  m_docaErr_prop = docaErrProp;
  m_poca = poca;
  m_pocaVector = pocaVector;
}  

/// set all the values
void AcdPocaDataV2::copy(const AcdPocaDataV2& other)                            
{
  set(other.getVolume(),other.getRegion(),other.getArcLength(),
      other.getDoca(),other.getDocaErrProj(),other.getDocaErrProp(),
      other.getPoca(),other.getPocaVector());         
}  


/// reset all the values to their default
void AcdPocaDataV2::clear() 
{
  m_volume = -1;
  m_region = -1;
  m_arcLength = 0.;
  m_doca = 0.;
  m_docaErr_proj = 0.;
  m_docaErr_prop = 0.;
  m_poca.SetXYZ(0.,0.,0.);
  m_pocaVector.SetXYZ(0.,0.,0.);
}

/// Print out this structure
void AcdPocaDataV2::print(Option_t *option) const 
{
  using namespace std;    
  cout << "AcdPocaDataV2 " << m_volume << ':' << m_region 
       << "s: " << m_arcLength << "  d: " << m_doca << " +- " << m_docaErr_proj << " +- " <<  m_docaErr_prop
       << "  p: [" << m_poca.X() << ' ' << m_poca.Y() << ' ' << m_poca.Z() << "] "
       << "  v: (" << m_pocaVector.X() << ' ' << m_pocaVector.Y() << ' ' << m_pocaVector.Z() << ") " << endl;
}

void AcdPocaDataV2::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  return;
}

Bool_t AcdPocaDataV2::CompareInRange( const AcdPocaDataV2& /* other */, const std::string & /* name */ ) const {
  return kTRUE;
}

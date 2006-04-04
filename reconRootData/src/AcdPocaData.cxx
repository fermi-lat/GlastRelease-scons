#include "reconRootData/AcdPocaData.h"
#include "Riostream.h"

ClassImp(AcdPocaData) ;

/// Default constructor.  Set everything to default values
AcdPocaData::AcdPocaData() {
  clear();
}

/// Constructor for use in transient -> persistent conversion 
/// Takes arguements as they are stored in ROOT
AcdPocaData::AcdPocaData(Float_t arcLength, Float_t doca, Float_t docaErr, 
			 const TVector3& poca, const TVector3& pocaVector){
  //Clear();
  set(arcLength,doca,docaErr,poca,pocaVector);
}

/// Copy constructor
AcdPocaData::AcdPocaData(const AcdPocaData& other)
{
  //Clear();
  set(other.getArcLength(),other.getDoca(),other.getDocaErr(),
      other.getPoca(),other.getPocaVector());	     
}

/// Assignment operator
AcdPocaData& AcdPocaData::operator=(const AcdPocaData& other)
{
  if ( this == &other ) return *this;
  set(other.getArcLength(),
      other.getDoca(),other.getDocaErr(),
      other.getPoca(),other.getPocaVector());
  return *this;
}


/// set all the values
void AcdPocaData::set(Float_t arcLength, Float_t doca, Float_t docaErr,
		      const TVector3& poca, const TVector3& pocaVector)			    
{
  m_arcLength = arcLength;
  m_doca = doca;
  m_docaErr = docaErr;
  m_poca = poca;
  m_pocaVector = pocaVector;
}  

/// set all the values
void AcdPocaData::copy(const AcdPocaData& other)			    
{
  m_arcLength = other.m_arcLength;
  m_doca = other.m_doca;
  m_docaErr = other.m_docaErr;
  m_poca = other.m_poca;
  m_pocaVector = other.m_pocaVector;
}  


/// reset all the values to their default
void AcdPocaData::clear() 
{
  m_arcLength = 0.;
  m_doca = 0.;
  m_docaErr = 0.;
  m_poca.SetXYZ(0.,0.,0.);
  m_pocaVector.SetXYZ(0.,0.,0.);
}

/// Print out this structure
void AcdPocaData::print(Option_t *option) const 
{
  using namespace std;    
  cout << "s: " << m_arcLength << "  d: " << m_doca << " +- " << m_docaErr 
       << "  p: [" << m_poca.X() << ' ' << m_poca.Y() << ' ' << m_poca.Z() << "] "
       << "  v: (" << m_pocaVector.X() << ' ' << m_pocaVector.Y() << ' ' << m_pocaVector.Z() << ") " << endl;
}

void AcdPocaData::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ ) {
  return;
}

Bool_t AcdPocaData::CompareInRange( const AcdPocaData& /* other */, const std::string & /* name */ ) const {
  return kTRUE;
}

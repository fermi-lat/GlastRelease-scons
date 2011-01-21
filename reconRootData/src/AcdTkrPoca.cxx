#include "reconRootData/AcdTkrPoca.h"
#include "Riostream.h"

ClassImp(AcdTkrPoca) ;

/// Default constructor.  Set everything to default values
AcdTkrPoca::AcdTkrPoca() 
  :TObject(){
  Clear("");
}

/// Constructor for use in transient -> persistent conversion 
/// Takes arguements as they are stored in ROOT
AcdTkrPoca::AcdTkrPoca(const AcdId& acdId, Int_t tkIndex,
                       Float_t doca, Float_t docaErr, UInt_t docaRegion,
                       const TVector3& poca, const TkrTrackParams& paramsAtPoca)
  :TObject()
{
  //Clear();
  initialize(acdId,tkIndex,doca,docaErr,docaRegion,poca,paramsAtPoca);
}

/// Copy constructor
AcdTkrPoca::AcdTkrPoca(const AcdTkrPoca& other)
  :TObject(other)
{
  //Clear();
  initialize(other.getId(),other.getTkrIndex(),
             other.getDoca(),other.getDocaErr(),other.getDocaRegion(),
             other.getPoca(),other.getParamsAtPoca());
}

/// Assignment operator
AcdTkrPoca& AcdTkrPoca::operator=(const AcdTkrPoca& other)
{
  if ( this == &other ) return *this;
  initialize(other.getId(),other.getTkrIndex(),
             other.getDoca(),other.getDocaErr(),other.getDocaRegion(),
             other.getPoca(),other.getParamsAtPoca());
  return *this;
}


/// set all the values
void AcdTkrPoca::initialize(const AcdId& acdId, Int_t tkIndex, 
                            Float_t doca, Float_t docaErr, UInt_t docaRegion,
                            const TVector3& poca, const TkrTrackParams& paramsAtPoca) {  
  m_id = acdId;
  m_tkIndex = tkIndex;
  m_doca = doca;
  m_docaErr = docaErr;
  m_docaRegion = docaRegion;
  m_poca = poca;
  m_paramsAtPoca = paramsAtPoca;
}  


/// reset all the values to their default
void AcdTkrPoca::Clear(Option_t* /* option */) 
{
  m_id = AcdId();
  m_tkIndex = -1;
  m_doca = 0.;
  m_docaErr = 0.;
  m_docaRegion = NONE_TILE;
  m_poca.SetXYZ(0.,0.,0.);
  m_paramsAtPoca = TkrTrackParams();
}

/// Print out this structure
void AcdTkrPoca::Print(Option_t *option) const 
{
  TObject::Print(option);
  using namespace std;    
}


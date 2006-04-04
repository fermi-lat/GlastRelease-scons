#include "reconRootData/AcdTkrLocalCoords.h"
#include "Riostream.h"

#include "TMatrixD.h"

ClassImp(AcdTkrLocalCoords)

AcdTkrLocalCoords::AcdTkrLocalCoords() 
{
  clear();
}


AcdTkrLocalCoords::AcdTkrLocalCoords(const float position[2], float pathLength, float cosTheta, 
				     int region, const TMatrixD& covMatrix)
{
  set(position,pathLength,cosTheta,region,covMatrix);
}

AcdTkrLocalCoords::AcdTkrLocalCoords(const AcdTkrLocalCoords& other)
{
  copy(other);
}

void AcdTkrLocalCoords::clear()
{
  
  m_activeX = 0.;
  m_activeY = 0.;    
  m_pathLength = 0.;
  m_cosTheta = 0.;
  m_XXCov = 0.;
  m_YYCov = 0.;
  m_XYCov = 0.;  

  m_region = 0;
}

void AcdTkrLocalCoords::print(Option_t * /*option*/ ) const 
{
  using namespace std;
  cout << "Active:  {" << m_activeX << ' ' << m_activeY << "}  path: "
       << m_pathLength << "  angle: " << m_cosTheta 
       << "  cov: <" << m_XXCov << ' ' << m_YYCov << ' ' << m_XYCov << "> " 
       << (int)m_region << endl;
}

AcdTkrLocalCoords& AcdTkrLocalCoords::operator=(const AcdTkrLocalCoords& params)
{
  copy(params);
  return *this;
}

void AcdTkrLocalCoords::set(const float position[2], float pathLength, float cosTheta, 
			    int region, const TMatrixD& CovMatrix)
{
  m_activeX = position[0];
  m_activeY = position[1];
  m_XXCov = CovMatrix[0][0];
  m_YYCov = CovMatrix[1][1];
  m_XYCov = CovMatrix[0][1];
  m_pathLength = pathLength;
  m_cosTheta = cosTheta;
  m_region = region;
}  

void AcdTkrLocalCoords::copy(const AcdTkrLocalCoords& other)
{
  m_activeX = other.getActiveX();
  m_activeY = other.getActiveY();
  m_XXCov = other.getLocalXXCov();
  m_YYCov = other.getLocalYYCov();
  m_XYCov = other.getLocalXYCov();
  m_pathLength = other.getPathLength();
  m_cosTheta = other.getCosTheta();
  m_region = other.getRegion();
}  


//======================================================
// For Unit Tests
//======================================================

void AcdTkrLocalCoords::Fake( Int_t /* ievent */, UInt_t /*rank*/, Float_t /*randNum*/ ) {

   // NOT TESTED YET

 }

Bool_t AcdTkrLocalCoords::CompareInRange( const AcdTkrLocalCoords & /*ref*/, const std::string & name ) const {

    bool result = true ;

    // NOT TESTED YET
    //result = ... && result ;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<Class_Name()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}


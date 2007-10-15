#include "reconRootData/AcdTkrIntersection.h"
#include "Riostream.h"
#include <cmath>

#include "TVector3.h"

ClassImp(AcdTkrIntersection)

AcdTkrIntersection::AcdTkrIntersection() 
  :TObject(){
    Clear();
}


AcdTkrIntersection::AcdTkrIntersection(const AcdId& tileId, Int_t trackIndex, 
				       const TVector3& globalPosition, 
				       const Double_t localPosition[2], const TMatrixD& localCovMatrix,
				       Double_t arcLengthToIntersection, Double_t pathLengthInTile, 
				       UChar_t tileHit, 
				       Double_t cosTheta)
  :TObject()
{
    initialize(tileId,trackIndex,
	       globalPosition, 
	       localPosition,localCovMatrix,
	       arcLengthToIntersection,pathLengthInTile,tileHit,cosTheta);
}

AcdTkrIntersection::AcdTkrIntersection(const AcdTkrIntersection& params)
  :TObject(params)
{
    initialize(params);
}

void AcdTkrIntersection::Clear(Option_t* /* option */) 
{
  
  m_tileId.Clear(); 
  m_location.Clear();

  m_trackIndex = -1;
  m_localX = 0.;
  m_localY = 0.;    
  m_localXXCov = 0.;
  m_localYYCov = 0.;
  m_localXYCov = 0.;  
  m_arcLengthToIntersection = 0.;  
  m_pathlengthInTile = 0.;

  m_tileHit = 0;
}

void AcdTkrIntersection::Print(Option_t *option) const 
{
    using namespace std;
    double localXErr = m_localXXCov > 0. ? sqrt(m_localXXCov) : 0;
    double localYErr = m_localYYCov > 0. ? sqrt(m_localYYCov) : 0.;
    double correl =  m_localXYCov / ( localXErr * localYErr );          
    std::cout << "AcdTkrIntersection.  id: " << m_tileId.getId() << "  s: " << m_arcLengthToIntersection 
	      << "  x: [" << m_location.X() << ',' << m_location.Y() << ',' <<  m_location.Z()
	      << "]  local: {" << m_localX << ',' << m_localY  
	      << "}  cov: <" << localXErr << ',' << localYErr << ',' << correl 
	      << ">  path: "<< m_pathlengthInTile << "  angle: " << m_cosTheta 
	      << "  hit: "
	      << (int)(m_tileHit)
	      << std::endl;
}

AcdTkrIntersection& AcdTkrIntersection::operator=(const AcdTkrIntersection& params)
{
    initialize(params);

    return *this;
}

void AcdTkrIntersection::initialize(const AcdId& tileId, Int_t trackIndex,
				    const TVector3& globalPostion, 
				    const Double_t localPosition[2], const TMatrixD& localCovMatrix,
				    Double_t arcLengthToIntersection, Double_t pathLengthInTile, 
				    UChar_t tileHit, Double_t cosTheta) {
  
  m_tileId = tileId;
  m_trackIndex = trackIndex;

  m_location = globalPostion;  
  m_localX = localPosition[0];
  m_localY = localPosition[1];
  m_localXXCov = localCovMatrix[0][0];
  m_localYYCov = localCovMatrix[1][1];
  m_localXYCov = localCovMatrix[0][1];
  m_arcLengthToIntersection = arcLengthToIntersection;  
  m_pathlengthInTile = pathLengthInTile;
  
  m_tileHit = tileHit;
  m_cosTheta = cosTheta;
}  

void AcdTkrIntersection::initialize(const AcdTkrIntersection& params)
{
  const Double_t local[2] = { params.getLocalX(), params.getLocalY() };
  TMatrixD localCov(2,2);
  localCov[0][0] = params.getLocalXXCov();
  localCov[0][1] = params.getLocalXYCov();
  localCov[1][0] = params.getLocalXYCov();
  localCov[1][1] = params.getLocalYYCov();
  initialize(params.getTileId(), params.getTrackIndex(), 
	     params.getGlobalPosition(),
	     local,localCov,
	     params.getArcLengthToIntersection(),params.getPathLengthInTile(),params.tileHit(),params.getCosTheta());
}  


//======================================================
// For Unit Tests
//======================================================

void AcdTkrIntersection::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

   // NOT TESTED YET

 }

Bool_t AcdTkrIntersection::CompareInRange( const AcdTkrIntersection & /*ref*/, const std::string & name ) const {

    bool result = true ;

    // NOT TESTED YET
    //result = ... && result ;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}


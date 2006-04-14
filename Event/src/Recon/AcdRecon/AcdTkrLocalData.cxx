// File and Version information:
// $Header$
//
//  Implementation file of AcdTkrPoca and AcdTkrPocaCol classes
//  
// Authors:
//
//    Eric Charles
//
//

#include "Event/Recon/AcdRecon/AcdTkrLocalData.h"

#include "GaudiKernel/MsgStream.h"
#include "CLHEP/Matrix/Matrix.h"


using namespace Event;

AcdTkrLocalCoords::AcdTkrLocalCoords() {
  ini();
}

AcdTkrLocalCoords::AcdTkrLocalCoords(const float localPosition[2], float pathLength, float cosTheta, 
				     int region, const CLHEP::HepMatrix& localCovMatrix) {
  set(localPosition,pathLength,cosTheta,region,localCovMatrix);
}

AcdTkrLocalCoords::AcdTkrLocalCoords(const AcdTkrLocalCoords& other) {
  CLHEP::HepMatrix localCov(2,2);
  localCov[0][0] = m_localXXCov;
  localCov[1][1] = m_localYYCov;
  localCov[0][1] = localCov[1][0] = m_localXYCov;
  float position[2];
  position[0] = other.m_activeX; position[0] = other.m_activeY;
  set(position,other.m_pathLength,other.m_cosTheta,other.m_region,localCov);
}

void AcdTkrLocalCoords::set(const float localPosition[2], float pathLength, float cosTheta, 
			    int region, const CLHEP::HepMatrix& localCov) {

  m_activeX = localPosition[0];
  m_activeY = localPosition[1];
  
  m_pathLength = pathLength;  
  m_cosTheta = cosTheta;  
  m_region = region;    
  
  m_localXXCov = localCov[0][0];
  m_localYYCov = localCov[1][1]; 
  m_localXYCov = localCov[1][0]; 
}


void AcdTkrLocalCoords::copy(const AcdTkrLocalCoords& other) {
  m_activeX = other.getActiveX();
  m_activeY = other.getActiveY();
  
  m_pathLength = other.getPathLength();
  m_cosTheta = other.getCosTheta();
  m_region = other.getRegion();    
  
  m_localXXCov = other.getLocalXXCov();
  m_localYYCov = other.getLocalYYCov();
  m_localXYCov = other.getLocalXYCov();
}


void AcdTkrLocalCoords::writeOut(MsgStream& stream ) const
// Purpose: provide ascii output of some data members for
//          debugging purposes
// Input:
//        stream - Gaudi message stream
{

  double localXErr = sqrt(m_localXXCov);
  double localYErr = sqrt(m_localYYCov);
  double correl =  m_localXYCov / ( localXErr * localYErr );
  stream << "Active: [" << m_activeX << ',' << m_activeY 
	 << "].  Cov: {" << localXErr << ',' << localYErr << ',' << correl 
	 << "}.  Path: " << m_pathLength 
	 << ".  Normal: " << m_cosTheta
	 << ".  Region: " << m_region
	 << ".  ";
}



void AcdTkrLocalCoords::ini()
// Purpose: reset all data members to 0
//
{
  m_activeX = 0.;
  m_activeY = 0.;
  
  m_pathLength = 0.;
  m_cosTheta = 0.;
  m_region = 0;
  
  m_localXXCov = 0.;
  m_localYYCov = 0.;
  m_localXYCov = 0.;
}

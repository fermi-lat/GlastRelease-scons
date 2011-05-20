// File and Version information:
// $Header$
//
//  Implementation file of AcdTkrIntersection and AcdTkrIntersectionCol classes
//  
// Authors:
//
//    Eric Charles
//
//
#include "Event/Recon/AcdRecon/AcdTkrIntersection.h"

#include "GaudiKernel/MsgStream.h"
#include "CLHEP/Matrix/Matrix.h"


using namespace Event;

AcdTkrIntersection::AcdTkrIntersection(){
  ini();
}

AcdTkrIntersection::AcdTkrIntersection(const idents::AcdId& acdId, int trackIndex,
                                       const Point& globalPosition, 
                                       const double localPosition[2], const CLHEP::HepMatrix& localCovMatrix,
                                       double arcLengthToIntersection, double pathLengthInTile,
                                       unsigned char tileHit, double cosTheta) {
  set(acdId,trackIndex, 
      globalPosition,localPosition, localCovMatrix,
      arcLengthToIntersection, pathLengthInTile,tileHit,cosTheta);
}

void AcdTkrIntersection::set(const idents::AcdId& acdId, int trackIndex, 
                             const Point& globalPosition, 
                             const double localPosition[2], const CLHEP::HepMatrix& localCovMatrix,
                             double arcLengthToIntersection, double pathLengthInTile,
                             unsigned char tileHit, double cosTheta) {
  m_tileId = acdId;
  m_trackIndex = trackIndex;
  m_location = globalPosition;
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
        

void AcdTkrIntersection::writeOut(MsgStream& stream) const

// Purpose: provide ascii output of some data members for
//          debugging purposes
// Input:
//        stream - Gaudi message stream
{

  double localXErr = sqrt(m_localXXCov);
  double localYErr = sqrt(m_localYYCov);
  double correl =  m_localXYCov / ( localXErr * localYErr );

  stream << MSG::DEBUG
         << "AcdTkrIntersection.  Tile: " << m_tileId.id()
         << ".  Track: " << m_trackIndex
         << ".  Global: (" << m_location.getX() << ',' << m_location.getY() << ',' << m_location.getZ()
         << ").  Local: [" << m_localX << ',' << m_localY
         << "].  Cov: {" << localXErr << ',' << localYErr << ',' << correl 
         << "}.  Arc: " << m_arcLengthToIntersection
         << ".  Path: " << m_pathlengthInTile
         << ".  HitMask: " << (int)m_tileHit
         << endreq;
}



void AcdTkrIntersection::ini()
// Purpose: reset all data members to 0
//
{
  
  m_tileId = idents::AcdId();
  m_trackIndex = -1;
  m_location.set(0.,0.,0.);
  m_localX = 0.;
  m_localY = 0.;
  
  m_localXXCov = 0.;
  m_localYYCov = 0.;
  m_localXYCov = 0.;    

  m_arcLengthToIntersection = 0.;
  m_pathlengthInTile = 0.;

  m_tileHit = 0;
  m_cosTheta = 0;
}

AcdTkrIntersectionCol::AcdTkrIntersectionCol(const std::vector<AcdTkrIntersection*>& acdTkrIntersections) {
  for ( std::vector<AcdTkrIntersection*>::const_iterator itr = acdTkrIntersections.begin();
        itr != acdTkrIntersections.end(); itr++ ) {
    AcdTkrIntersection* iSect = const_cast<AcdTkrIntersection*>(*itr);
    add(iSect);
  }
}

void AcdTkrIntersectionCol::delIntersections() 

//Purpose: delete all AcdTkrIntersection object from memory

{
    int nInter = num();
    for (int iIn = 0; iIn < nInter; iIn++) {
        delete operator[](iIn);
    }
    clear();
}

void AcdTkrIntersectionCol::ini()

//Purpose:  delete all pointers to clusters
// from collection 

{
    clear();
}


void AcdTkrIntersectionCol::writeOut(MsgStream& stream) const

// Purpose: provide symbolic output of some data members
//          of all clusters in collection for debugging purposes
//
// Input:
//        stream - Gaudi message stream
{
    
  // loop over all clusters
  for (unsigned int i = 0; i < size();i++) {
    
    // call the writeOut() method for each cluster
    (operator[](i))->writeOut(stream);
  }
  
}

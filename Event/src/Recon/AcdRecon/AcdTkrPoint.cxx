// File and Version information:
// $Header$
//
//  Implementation file of AcdTkrPoint and AcdTkrPointCol classes
//  
// Authors:
//
//    Eric Charles
//
//

#include "Event/Recon/AcdRecon/AcdTkrPoint.h"

#include "GaudiKernel/MsgStream.h"
#include "CLHEP/Matrix/Matrix.h"


using namespace Event;

AcdTkrPoint::AcdTkrPoint(){
  ini();
}

AcdTkrPoint::AcdTkrPoint(const AcdTkrPoint& other)
  :AcdTkrLocalCoords(other),
   m_trackIndex(other.m_trackIndex){
}

/// Assignment operator
Event::AcdTkrPoint& 
AcdTkrPoint::operator=(const Event::AcdTkrPoint& other) {
  m_trackIndex = other.getTrackIndex();
  AcdTkrLocalCoords::copy(other);
  return *this;
}

AcdTkrPoint::AcdTkrPoint( float arcLength, int volume,
                          const HepPoint3D& global, const Event::TkrTrackParams& /* params */ )
  :AcdTkrLocalCoords(volume,arcLength,global),
   m_trackIndex(-1){
}


AcdTkrPoint::AcdTkrPoint(int trackIndex,
                         int volumePlane, float arcLengthToPlane, float cosTheta, 
                         const HepPoint3D& global, const float localPosition[2], 
                         const CLHEP::HepSymMatrix& localCovProj, const CLHEP::HepSymMatrix& localCovProp)
  :AcdTkrLocalCoords(volumePlane,arcLengthToPlane,cosTheta,
                     global,localPosition,localPosition,
                     localCovProj,localCovProp),
   m_trackIndex(trackIndex){
}

void AcdTkrPoint::set(int trackIndex) {
  m_trackIndex = trackIndex;
}
        

void AcdTkrPoint::writeOut(MsgStream& stream) const

// Purpose: provide ascii output of some data members for
//          debugging purposes
// Input:
//        stream - Gaudi message stream
{

  stream << MSG::DEBUG
         << "AcdTkrPoint:  " << m_trackIndex << std::endl;
  AcdTkrLocalCoords::writeOut(stream);
  stream << endreq;
}



void AcdTkrPoint::ini()
// Purpose: reset all data members to 0
//
{  
  m_trackIndex = -1;
  AcdTkrLocalCoords::ini();
}

AcdTkrPointCol::AcdTkrPointCol(const std::vector<AcdTkrPoint*>& acdTkrPoints) {
  for ( std::vector<AcdTkrPoint*>::const_iterator itr = acdTkrPoints.begin();
        itr != acdTkrPoints.end(); itr++ ) {
    AcdTkrPoint* iSect = const_cast<AcdTkrPoint*>(*itr);
    add(iSect);
  }
}

void AcdTkrPointCol::delTkrPoints() 

//Purpose: delete all AcdTkrPoint object from memory

{
    int nInter = num();
    for (int iIn = 0; iIn < nInter; iIn++) {
        delete operator[](iIn);
    }
    clear();
}

void AcdTkrPointCol::ini()

//Purpose:  delete all pointers to clusters
// from collection 

{
    clear();
}


void AcdTkrPointCol::writeOut(MsgStream& stream) const

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

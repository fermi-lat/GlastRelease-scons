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

AcdTkrPoint::AcdTkrPoint(double arcLength, int face, 
			 const Point& point, const Event::TkrTrackParams& paramsAtPoint){
  set(arcLength,face,point,paramsAtPoint);
}

void AcdTkrPoint::set(double arcLength, int face, 
		      const Point& point, const Event::TkrTrackParams& paramsAtPoint) {
  m_arcLength = arcLength;
  m_face = face;
  m_point = point;
  m_paramsAtPoint = paramsAtPoint;
}
        

void AcdTkrPoint::writeOut(MsgStream& stream) const

// Purpose: provide ascii output of some data members for
//          debugging purposes
// Input:
//        stream - Gaudi message stream
{

  stream << MSG::DEBUG
	 << "AcdTkrPoint." 
	 << "  Arc: " << m_arcLength  
	 << ".  Face: " << m_face
	 << ".  Global: (" << m_point.x() << ',' << m_point.y() << ',' << m_point.z() << ")"
    /* << m_paramsAtPoint */
	 << endreq;
}



void AcdTkrPoint::ini()
// Purpose: reset all data members to 0
//
{  
  m_arcLength = 0.;
  m_face = -1;
  m_point = Point();
  m_paramsAtPoint = Event::TkrTrackParams();
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

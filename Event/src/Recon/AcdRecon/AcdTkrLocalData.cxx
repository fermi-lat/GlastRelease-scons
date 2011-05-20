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

#include "geometry/Point.h"


using namespace Event;

AcdTkrLocalCoords::AcdTkrLocalCoords():
  m_volume(-1),
  m_arcLengthPlane(0.),
  m_cosTheta(0),
  m_global(0.,0.,0.),
  m_localCovProj(2,0),
  m_localCovProp(2,0){
  m_local[0] = m_local[1] = 0.;
  m_active[0] = m_active[1] = 0.;
}

AcdTkrLocalCoords::AcdTkrLocalCoords(int volume, float arcLengthPlane, 
                                     const HepPoint3D& global )
  :m_volume(volume),
   m_arcLengthPlane(arcLengthPlane),
   m_cosTheta(0),
   m_global(global),
   m_localCovProj(2,0),
   m_localCovProp(2,0){
  m_local[0] = m_local[1] = 0.;
  m_active[0] = m_active[1] = 0.;
}
      
AcdTkrLocalCoords::AcdTkrLocalCoords(int volume, float arcLengthPlane, float cosTheta, 
                                     const HepPoint3D& global, 
                                     const float localPosition[2], const float active[2], 
                                     const CLHEP::HepSymMatrix& localCovProj, const CLHEP::HepSymMatrix& localCovProp)
  :m_volume(volume),
   m_arcLengthPlane(arcLengthPlane),
   m_cosTheta(cosTheta),
   m_global(global),
   m_localCovProj(localCovProj),
   m_localCovProp(localCovProp){
  m_local[0] = localPosition[0];
  m_local[1] = localPosition[1];  
  m_active[0] = active[0];
  m_active[1] = active[1];
}
    
AcdTkrLocalCoords::AcdTkrLocalCoords(float arcLength, float cosTheta, 
                                     const HepPoint3D& global, 
                                     const double localPosition[2], 
                                     const CLHEP::HepSymMatrix& planeError)
  :m_volume(-1),
   m_arcLengthPlane(arcLength),
   m_cosTheta(cosTheta),
   m_global(global),
   m_localCovProj(planeError),
   m_localCovProp(planeError){
  m_local[0] = localPosition[0];
  m_local[1] = localPosition[1];  
  m_active[0] = localPosition[0];
  m_active[1] = localPosition[1];
}
    


AcdTkrLocalCoords::AcdTkrLocalCoords(const AcdTkrLocalCoords& other)
  :m_volume(other.m_volume),
   m_arcLengthPlane(other.m_arcLengthPlane),
   m_cosTheta(other.m_cosTheta),
   m_global(other.m_global),
   m_localCovProj(other.m_localCovProj),
   m_localCovProp(other.m_localCovProp){
  m_local[0] = other.m_local[0];
  m_local[1] = other.m_local[1]; 
  m_active[0] = other.m_active[0];
  m_active[1] = other.m_active[1];
}

void AcdTkrLocalCoords::setLocalData(int volume, float arcLengthPlane, float cosTheta, 
                                     const HepPoint3D& global, const float localPosition[2], const float active[2],
                                     const CLHEP::HepSymMatrix& localCovProj, const CLHEP::HepSymMatrix& localCovProp) {
  
  m_volume =  volume;
  m_arcLengthPlane = arcLengthPlane;
  m_cosTheta = cosTheta;
  m_global = global;
  m_local[0] = localPosition[0];
  m_local[1] = localPosition[1];
  m_active[0] = active[0];
  m_active[1] = active[1];
  m_localCovProj = localCovProj;
  m_localCovProp = localCovProp;
}

/// set everything at once, old version
void AcdTkrLocalCoords::setLocalData(const float localPosition[2],  
                                     float /* pathLength */, float cosTheta, 
                                     int region, const CLHEP::HepSymMatrix& planeError) {
  m_volume = region;
  m_local[0] = localPosition[0];
  m_local[1] = localPosition[1];
  m_active[0] = localPosition[0];
  m_active[1] = localPosition[1];
  m_cosTheta = cosTheta;
  m_localCovProj = planeError;
  m_localCovProp = planeError;
}

/// set stuff from old version of AcdTkrPoint
void AcdTkrLocalCoords::setLocalData(float arcLength, int face, 
                                     const Point& point, const Event::TkrTrackParams& /*params*/) {
  m_volume = face;
  m_arcLengthPlane = arcLength;
  m_global.set( point.x(), point.y(), point.z());
}


void AcdTkrLocalCoords::copy(const AcdTkrLocalCoords& other) {
  m_volume = other.getLocalVolume();
  m_arcLengthPlane = other.getArclengthToPlane();
  m_cosTheta = other.getCosTheta();
  m_global = other.getGlobalPosition();
  m_local[0] = other.getLocalX();
  m_local[1] = other.getLocalY();
  m_active[0] = other.getActiveX();
  m_active[1] = other.getActiveY();
  m_localCovProj = other.getLocalCovProj();
  m_localCovProp = other.getLocalCovProp();
}


void AcdTkrLocalCoords::writeOut(MsgStream& stream ) const
// Purpose: provide ascii output of some data members for
//          debugging purposes
// Input:
//        stream - Gaudi message stream
{

  stream << "  LocalVol: " << m_volume 
         << "  Global. [" << m_global.x() << ',' << m_global.y() << ',' << m_global.z() 
         << "]. Local: [" << m_local[0] << ',' << m_local[1] 
         << "].  S: " << m_arcLengthPlane
         << ".   Angle: " << m_cosTheta
         << ".   CovProjected:  {" << m_localCovProj
         << "}.   CovPropagated:  {" << m_localCovProp
         << "}.";
}



void AcdTkrLocalCoords::ini()
// Purpose: reset all data members to 0
//
{
  m_arcLengthPlane = 0.;
  m_cosTheta = 0.;
  m_global.set(0.,0.,0.);
  m_local[0] = 0.;
  m_local[1] = 0.;
  m_active[0] = 0.;
  m_active[1] = 0.;
  m_localCovProj = CLHEP::HepSymMatrix(2,0);
  m_localCovProp = CLHEP::HepSymMatrix(2,0);  
}

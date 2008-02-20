// File and Version information:
// $Header$
//
//  Implementation file of AcdRibbonDim 
//  
// Authors:
//
//    Eric Charles
//
//

#include "AcdUtil/AcdRibbonDim.h"

#include "CLHEP/Geometry/Transform3D.h"

/// Constructor: takes the ribbon id, the volume id, and the detector service
AcdRibbonDim::AcdRibbonDim(const idents::AcdId& acdId, IAcdGeometrySvc& acdGeomSvc) 
  :m_acdId(acdId),
   m_acdGeomSvc(acdGeomSvc){
  m_sc = getVals();
}  

/// this function access the detector service to get the geometry information
StatusCode AcdRibbonDim::getVals() {

  StatusCode sc = StatusCode::SUCCESS;

  // First we load up the rays
  m_minusSideRays.clear();
  m_topRays.clear();
  m_plusSideRays.clear();
  m_halfWidth = m_acdGeomSvc.ribbonHalfWidth();
  bool isOk = m_acdGeomSvc.fillRibbonData(m_acdId,m_minusSideRays,m_topRays,m_plusSideRays,
					  m_minusSideTransform,m_topTransform,m_plusSideTransform);
  if ( ! isOk ) {
    return StatusCode::FAILURE;
  }
  double fullLength = calcRibbonLength();
  if ( fullLength <= 0. ) {
    return StatusCode::FAILURE;
      }
  m_halfLength = fullLength/2.;
  return sc;
}


void AcdRibbonDim::toLocal(const HepPoint3D& global, int isegment, HepPoint3D& local) const {
  switch (isegment) {
  case 0: 
    local = m_minusSideTransform * global;
    break;
  case 1: 
    local = m_topTransform * global;
    break;
  case 2: 
    local = m_plusSideTransform * global;
    break;    
  }
}

bool AcdRibbonDim::setEdgeRay(int iSeg, HepPoint3D& start, HepVector3D& vector) const {
  const Ray* firstRay(0);
  const Ray* lastRay(0);
  switch (iSeg) {
  case 0: 
    firstRay = &(m_minusSideRays[0]);
    lastRay = &(m_minusSideRays.back());
    break;
  case 1: 
    firstRay =  &(m_topRays[0]);
    lastRay = &(m_topRays.back());
    break;
  case 2: 
    firstRay = &(m_plusSideRays[0]);
    lastRay = &(m_plusSideRays.back());
    break;    
  }
  const Point& startPos = firstRay->position();
  Point endPos = lastRay->position(lastRay->getArcLength());
  start.set(startPos.x(),startPos.y(),startPos.z());
  switch (iSeg) {
  case 0: 
  case 2:
    // up the sides.  only z counts
    vector.set(0.,0.,endPos.z()-startPos.z());
    break;
  case 1: 
    // across the top. only x or y count
    vector.set(endPos.x()-startPos.x(),endPos.y()-startPos.y(),0);
    break;    
  }
  return true;
}

const HepTransform3D& AcdRibbonDim::transform(int iVol) const {
  switch ( iVol ) {
  case 0:
    return m_minusSideTransform;
  case 1:
    return m_topTransform;
  case 2:
    return m_plusSideTransform;
  }
  static const HepTransform3D nullTrans;
  return nullTrans;
}


double AcdRibbonDim::calcRibbonLength() {
   double len(0);
   std::vector<Ray>::const_iterator itr = m_minusSideRays.begin();
   
   for ( itr = m_minusSideRays.begin(); itr != m_minusSideRays.end(); itr++ ) {
     len += itr->getArcLength();
   }
   for ( itr = m_topRays.begin(); itr != m_topRays.end(); itr++ ) {
     len += itr->getArcLength();
   }
   for ( itr = m_plusSideRays.begin(); itr != m_plusSideRays.end(); itr++ ) {
     len += itr->getArcLength();
   }
   return len;  
}


bool AcdRibbonDim::getRibbonLengthAndBin(const idents::VolumeIdentifier &volId, const HepPoint3D& point,
					 double& /* ribbonLength */, int& ribbonBin) const {
  
  // sanity check
  idents::AcdId checkId(volId);
  if ( m_acdId != checkId ) { 
    return false;
  }

  int orient = m_acdId.ribbonOrientation();
  int face = volId[1];  
  int segment = volId[5];
  int code = (1000*orient) + (100* face) + segment;

  switch ( code ) {
  case 5001: // top
    ribbonBin = 3;
    break;
  case 5101: // -x side, top
  case 5106: // -x side, upper short
    ribbonBin = 1;
    break;
  case 5102: // -x side, middle
  case 5107: // -x side, lower short
  case 5103: // -x side, bottom
    ribbonBin = 0;
    break;
  case 5301: // +x side, top
  case 5306: // +x side, upper short
    ribbonBin = 5;
    break;
  case 5302: // +x side, middle
  case 5303: // +x side, bottom
  case 5307: // +x side, lower short
    ribbonBin = 6;
    break;
  case 6001: // top -y side
  case 6008: // top -y side short
    ribbonBin = 2;
    break;
  case 6002: // top -y middle
    ribbonBin = point.y() < 0 ? 2 : 3;
    break;
  case 6009: // top -y middle short
  case 6003: // top center
  case 6010: // top +y middle short
    ribbonBin = 3;
    break;
  case 6004: // top +y middle
    ribbonBin = point.y() < 0 ? 3 : 4;
    break;
  case 6011: // top +y side short
  case 6005: // top +y side
    ribbonBin = 4;
    break;
  case 6209: // -y side, bend short
  case 6205: // -y side, upper short
  case 6201: // -y side, top
    ribbonBin = 1;
    break;
  case 6206: // -y side, middle short
  case 6202: // -y side, middle
  case 6207: // -y side, lower short
  case 6203: // -y side, bottom
    ribbonBin = 0;
    break;
  case 6409: // +y side, bend short
  case 6405: // +y side, upper short    
  case 6401: // +y side, top
    ribbonBin = 5;
    break;
  case 6406: // +y side, middle short
  case 6402: // +y side, middle
  case 6407: // +y side, lower short
  case 6403: // +y side, bottom
    ribbonBin = 6;
    break;
  default:
    return false;
  }
  return true;  
}

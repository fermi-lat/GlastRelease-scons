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
  static const int ribbonX = 5 /*, ribbonY = 6; */;

  // First we load up the rays
  m_minusSideRays.clear();
  m_topRays.clear();
  m_plusSideRays.clear();
  m_halfWidth = m_acdGeomSvc.ribbonHalfWidth();
  bool isOk = m_acdGeomSvc.fillRibbonRays(m_acdId,m_minusSideRays,m_topRays,m_plusSideRays);
  if ( ! isOk ) {
    return StatusCode::FAILURE;
  }
  
  int minusFace = m_acdId.ribbonOrientation() == ribbonX ? 1 : 2;
  int plusFace = m_acdId.ribbonOrientation() == ribbonX ? 3 : 4;
  
  isOk = m_acdGeomSvc.fillRibbonTransform(minusFace,m_minusSideRays[0],m_minusSideTransform);
  if ( ! isOk ) {
    return StatusCode::FAILURE;
  }
  isOk = m_acdGeomSvc.fillRibbonTransform(plusFace,m_plusSideRays[0],m_plusSideTransform);
  if ( ! isOk ) {
    return StatusCode::FAILURE;
  }
  isOk = m_acdGeomSvc.fillRibbonTransform(0,m_topRays[0],m_topTransform);
  if ( ! isOk ) {
    return StatusCode::FAILURE;
  }
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

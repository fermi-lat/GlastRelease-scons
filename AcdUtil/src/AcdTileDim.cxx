// File and Version information:
// $Header$
//
//  Implementation file of AcdTileDim 
//  
// Authors:
//
//    Eric Charles
//
//

#include "AcdUtil/AcdTileDim.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "CLHEP/Geometry/Transform3D.h"

/// Constructor: takes the tile id, the volume id, and the detector service
AcdTileDim::AcdTileDim(const idents::AcdId& acdId, const idents::VolumeIdentifier& volId, 
		       IGlastDetSvc &detSvc) 
  :m_acdId(acdId),
   m_nVol(1),
   m_detSvc(detSvc) {
  m_volId[0] = volId;
  m_volId[1] = idents::VolumeIdentifier();
  m_sc = getVals();
}  

AcdTileDim::AcdTileDim(const idents::AcdId& acdId, 
		       const idents::VolumeIdentifier& volIdMain, 
		       const idents::VolumeIdentifier& volIdOther, IGlastDetSvc &detSvc)
  :m_acdId(acdId),
   m_nVol(2),
   m_detSvc(detSvc){
  m_volId[0] = volIdMain;
  m_volId[1] = volIdOther;
  m_sc = getVals();
}

/// this function access the detector service to get the geometry information
StatusCode AcdTileDim::getVals() {
  
  /// get the tile dimensions
  std::string str;
  StatusCode sc = StatusCode::SUCCESS;

  m_shared[0] = -1;
  m_shared[1] = -1;

  m_sharedWidth[0] = 0.;
  m_sharedWidth[1] = 0.;

  for ( int iVol(0); iVol < m_nVol; iVol++ ) {

    const idents::VolumeIdentifier& volId = m_volId[iVol];
    sc = m_detSvc.getShapeByID(volId, &str, &(m_dim[iVol]));
    if ( sc.isFailure() ) {
      //log << MSG::WARNING << "Failed to retrieve Shape by Id" << endreq;
      return sc;
    } 

    /// get the tile position
    sc = m_detSvc.getTransform3DByID(volId, &(m_transform[iVol]));
    if (sc.isFailure() ) {
      //log << MSG::WARNING << "Failed to get trasnformation" << endreq;
      return sc;
    } 
    HepPoint3D center(0., 0., 0.);
    m_tileCenter[iVol] = m_transform[iVol] * center;

    /// calculate the corners of the tile
    sc = getCorners(m_dim[iVol],m_tileCenter[iVol],&(m_corners[iVol][0]));
    if (sc.isFailure() ) {
      //log << MSG::WARNING << "Failed to get trasnformation" << endreq;
      return sc;
    } 
  }

  if ( m_nVol == 2 ) {
    switch (m_acdId.id() ) {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
      m_shared[0] = 3;
      m_shared[1] = 1;
      m_sharedWidth[0] = m_dim[1][2];
      m_sharedWidth[1] = m_dim[0][1];
      break;
    case 40:
    case 41:
    case 42:
    case 43:
    case 44:
      m_shared[0] = 1;
      m_shared[1] = 1;
      m_sharedWidth[0] = m_dim[1][2];
      m_sharedWidth[1] = m_dim[0][1];
      break;
    default:
      return StatusCode::FAILURE;      
    }
  }

  return sc;
}


StatusCode AcdTileDim::getCorners(const std::vector<double> &dim, const HepPoint3D &center, HepPoint3D *corner) {

  // fill corner with the 4 corner points of this tile
  StatusCode sc = StatusCode::SUCCESS;
  
  unsigned int iCorner;
  // Ignore short dimension - only interested in 4 corners
  
  if ((dim[0] < dim[1]) && (dim[0] < dim[2])) { //X smallest - X Face
    for(iCorner = 0; iCorner<4; iCorner++) {
      corner[iCorner].setX(center.x());
      corner[iCorner].setY(center.y() + 
			   ( (iCorner < 2) ? -1 : 1) * dim[1]*0.5);
      corner[iCorner].setZ(center.z() +
			   ( ((iCorner == 0) || (iCorner==3)) ? -1 : 1) * dim[2]*0.5);
    }    
  } else if ((dim[1] < dim[0]) && (dim[1] < dim[2])) {  //Y smallest - Y FACE
    for(iCorner = 0; iCorner<4; iCorner++) {
      corner[iCorner].setY(center.y());
      corner[iCorner].setX(center.x() + 
			   ( (iCorner < 2) ? -1 : 1) * dim[0]*0.5);
      corner[iCorner].setZ(center.z() +
			   ( ((iCorner == 0) || (iCorner==3)) ? -1 : 1) * dim[2]*0.5);
    }    
  } else {
    for(iCorner = 0; iCorner<4; iCorner++) { //Z smallest - TOP
      corner[iCorner].setZ(center.z());
      corner[iCorner].setX(center.x() + 
			   ( (iCorner < 2) ? -1 : 1) * dim[0]*0.5);
      corner[iCorner].setY(center.y() +
			   ( ((iCorner == 0) || (iCorner==3)) ? -1 : 1) * dim[1]*0.5);
    }
  }  

  return sc;
} 

void AcdTileDim::toLocal(const HepPoint3D& global, HepPoint3D& local, int idx) {
  assert(idx < m_nVol);
  local = m_transform[idx] * global;
}

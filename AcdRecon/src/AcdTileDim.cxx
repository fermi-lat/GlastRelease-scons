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

//#include "AcdUtil/AcdTileDim.h"
#include "./AcdTileDim.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "CLHEP/Geometry/Transform3D.h"

/// Constructor: takes the tile id, the volume id, and the detector service
AcdTileDim::AcdTileDim(const idents::AcdId& acdId, const idents::VolumeIdentifier& volId, 
		       IGlastDetSvc &detSvc) 
  :m_acdId(acdId),
   m_volId(volId),
   m_detSvc(detSvc) {
  m_sc = getVals();
}  

/// this function access the detector service to get the geometry information
StatusCode AcdTileDim::getVals() {
  
  /// get the tile dimensions
  std::string str;
  StatusCode sc = StatusCode::SUCCESS;
  sc = m_detSvc.getShapeByID(m_volId, &str, &m_dim);
  if ( sc.isFailure() ) {
    //log << MSG::WARNING << "Failed to retrieve Shape by Id" << endreq;
    return sc;
  } 

  /// get the tile position
  HepTransform3D transform;
  sc = m_detSvc.getTransform3DByID(m_volId, &transform);
  if (sc.isFailure() ) {
    //log << MSG::WARNING << "Failed to get trasnformation" << endreq;
    return sc;
  } 
  HepPoint3D center(0., 0., 0.);
  m_tileCenter = transform * center;

  /// calculate the corners of the tile
  sc = getCorners(m_dim,m_tileCenter,m_corners);
  if (sc.isFailure() ) {
    //log << MSG::WARNING << "Failed to get trasnformation" << endreq;
    return sc;
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

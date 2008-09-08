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
AcdTileDim::AcdTileDim(const idents::AcdId& acdId,  
		       IAcdGeometrySvc& acdGeomSvc)
  :m_acdId(acdId),
   m_acdGeomSvc(acdGeomSvc){
  m_sc = getVals();
}  


/// this function access the detector service to get the geometry information
StatusCode AcdTileDim::getVals() {
  
  std::map<idents::AcdId, int>::const_iterator itrFind = m_acdGeomSvc.getAcdIdVolCountCol().find(m_acdId);
  if ( itrFind == m_acdGeomSvc.getAcdIdVolCountCol().end() ) {
    return StatusCode::FAILURE;
  }
  m_nVol = itrFind->second;

  bool isOk = m_acdGeomSvc.fillScrewHoleData(m_acdId,m_screwHoles);
  if ( ! isOk ) return StatusCode::FAILURE;

  for ( int iVol(0); iVol < m_nVol; iVol++ ) {
 
    m_dim[iVol].clear();
    isOk = m_acdGeomSvc.fillTileData(m_acdId,iVol,m_transform[iVol],m_dim[iVol],m_tileCenter[iVol],m_corners[iVol]);
    if ( ! isOk ) return StatusCode::FAILURE;

    HepTransform3D inverseT = m_transform[iVol].inverse();
    HepVector3D xV(1.,0.,0.); HepVector3D yV(0.,1.,0.);
    xV = inverseT * xV; yV = inverseT * yV;
    m_localFrameVectors[iVol] = HepMatrix(2,3);
    m_localFrameVectors[iVol](1,1) = xV.x(); m_localFrameVectors[iVol](1,2) = xV.y(); m_localFrameVectors[iVol](1,3) = xV.z(); 
    m_localFrameVectors[iVol](2,1) = yV.x(); m_localFrameVectors[iVol](2,2) = yV.y(); m_localFrameVectors[iVol](2,3) = yV.z(); 
  }
  isOk = m_acdGeomSvc.fillTileSharedEdgeData(m_acdId,m_dim[0],m_dim[1],
					     m_shared[0],m_shared[1],m_sharedWidth[0],m_sharedWidth[1]);
  if ( ! isOk ) return StatusCode::FAILURE;

  if ( m_dim[0].size() == 5 ) {
    m_trapezoid = true;
    m_xmean = m_dim[0][0] + m_dim[0][3];
    m_xmean /= 2.;
    m_xdiff = m_dim[0][3] - m_dim[0][0];
    m_xdiff /= 2.;
  } else {
    m_trapezoid = false;
    m_xmean = 0.;
    m_xdiff = 0.;
  }

  return StatusCode::SUCCESS;
}

StatusCode AcdTileDim::toLocalCoords(const AcdTileDim& /* dim */, 
				     int /* region */, const double& /* activeX */, const double& /* activeY */,
				     double& localX, double& localY){
  localX = localY = 0.;
  return StatusCode::SUCCESS;
}


double AcdTileDim::widthAtLocalY(int iVol, double localY) const {
  if ( m_trapezoid ) {
    return m_xmean + m_xdiff * ( localY / m_dim[iVol][1] );
  } 
  return m_dim[iVol][0];
}


void AcdTileDim::activeDistance(const HepPoint3D& localPoint, int iVol, double& activeX, double& activeY) const {

  activeY = (dim(iVol)[1]/2.) - fabs(localPoint.y());
  if ( ! m_trapezoid ) {
    activeX = (dim(iVol)[0]/2.) - fabs(localPoint.x());
  } else {
    double correctionAtLocalY = ( m_dim[0][4] * (  localPoint.y() / m_dim[iVol][1] ) ) / 2.;
    double half_width = widthAtLocalY(iVol,localPoint.y())/2.;
    activeX = half_width - fabs(localPoint.x() + correctionAtLocalY); 
  }

  if ( nVol() == 2 ) {
    int face = acdId().face();
    if ( face == 0 ) {

      // first make sure that it actually hit the tile piece
      if ( activeY < 0. ) return;
      // Ok, it hit the tile.  Correct for the fact that it bends around
      if ( iVol == 0 ) {
	// top, check for hitting near bent tile
	if ( sharedEdge(0) == 1 &&  localPoint.y() > 0 ) {
	  // hit +Y side of TOP tile near +Y edge,
	  // add the extra size of the bent piece
	  activeY += fabs(sharedWidth(0));
	} else if (  sharedEdge(0) == 3 &&  localPoint.y() <  0 ) {
	  // hit -Y side of TOP tile near -Y edge,
	  // add the extra size of the bent piece
	  activeY += fabs(sharedWidth(0));
	}  
      } else if ( iVol == 1 ) {	
	// side, check for hitting near top tile
	if ( sharedEdge(1) == 1 && localPoint.y() > 0 ) {
	  // hit upper part of BENT piece on -Y side ( local Y goes UP ) 
	  // is a shared piece.  but this is a short side, so take the distance to the
	  // other side of this volume
	  activeY = (dim(iVol)[1]/2.) + localPoint.y();
	} else if ( sharedEdge(1) == 3 && localPoint.y() < 0 ) {
	  // hit upper part of BENT piece on +Y side ( local Y goes DOWN )
	  // is a shared piece.  but this is a short side, so take the distance to the
	  // other side of this volume
	  activeY = (dim(iVol)[1]/2.) - localPoint.y();
	}
      } 
    }
  }
}

void AcdTileDim::toLocal(const HepPoint3D& global, HepPoint3D& local, int idx) const {
  assert(idx < m_nVol);
  local = m_transform[idx] * global;
}




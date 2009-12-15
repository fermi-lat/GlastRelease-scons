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

AcdTileSection::AcdTileSection()
  :m_trapezoid(false),
   m_xmean(0.),
   m_xdiff(0.),
   m_shared(-1),
   m_sharedWidth(0.){
}


double AcdTileSection::widthAtLocalY(double localY) const {
  if ( m_trapezoid ) {
    return m_xmean + m_xdiff * ( localY / m_dim[1] );
  } 
  return m_dim[0];
}


void AcdTileSection::activeDistance(const HepPoint3D& localPoint, int iVol, double& activeX, double& activeY) const {
  
  activeX = (widthAtLocalY(localPoint.y())/2.) - fabs(localPoint.x() - m_xdiff); 
  activeY = (m_dim[1]/2.) - fabs(localPoint.y());

  // If no shared edges, we are done
  if ( m_shared < 0 ) return;
  
  // first make sure that it actually hit the tile piece
  if ( activeY < 0. ) return;
  // Ok, it hit the tile.  Correct for the fact that it bends around
  if ( iVol == 0 ) {
    // top, check for hitting near bent tile
    if ( m_shared == 1 &&  localPoint.y() > 0 ) {
      // hit +Y side of TOP tile near +Y edge,
      // add the extra size of the bent piece
      activeY += fabs(m_sharedWidth);
    } else if (  m_shared == 3 &&  localPoint.y() <  0 ) {
      // hit -Y side of TOP tile near -Y edge,
      // add the extra size of the bent piece
      activeY += fabs(m_sharedWidth);
    }  
  } else if ( iVol == 1 ) {	
    // side, check for hitting near top tile
    if ( m_shared == 1 && localPoint.y() > 0 ) {
      // hit upper part of BENT piece on -Y side ( local Y goes UP ) 
      // is a shared piece.  but this is a short side, so take the distance to the
      // other side of this volume
      activeY = (m_dim[1]/2.) + localPoint.y();
    } else if ( m_shared == 3 && localPoint.y() < 0 ) {
      // hit upper part of BENT piece on +Y side ( local Y goes DOWN )
      // is a shared piece.  but this is a short side, so take the distance to the
      // other side of this volume
      activeY = (m_dim[1]/2.) - localPoint.y();
    }
  }
}

/// Constructor: takes the tile id, the volume id, and the detector service
AcdTileDim::AcdTileDim(const idents::AcdId& acdId,  
		       IAcdGeometrySvc& acdGeomSvc)
  :m_acdId(acdId),
   m_acdGeomSvc(acdGeomSvc),
   m_doneGaps(false){
  m_minusXGap[0] = 0.;
  m_minusXGap[1] = 0.;
  m_plusXGap[0] = 0.;
  m_plusXGap[1] = 0.;
  m_sc = getVals();
}  


AcdTileDim::~AcdTileDim() {
  for ( std::vector<AcdTileSection*>::iterator itr = m_sections.begin(); itr != m_sections.end(); itr++ ) {
    AcdTileSection* sect = *itr;
    delete sect;
  }
}


/// this function access the detector service to get the geometry information
StatusCode AcdTileDim::getVals() {
  
  std::map<idents::AcdId, int>::const_iterator itrFind = m_acdGeomSvc.getAcdIdVolCountCol().find(m_acdId);
  if ( itrFind == m_acdGeomSvc.getAcdIdVolCountCol().end() ) {
    return StatusCode::FAILURE;
  }
  int nVol = itrFind->second;

  bool isOk = m_acdGeomSvc.fillScrewHoleData(m_acdId,m_screwHoles);
  if ( ! isOk ) return StatusCode::FAILURE;

  for ( int iVol(0); iVol < nVol; iVol++ ) {

    AcdTileSection* sect = new AcdTileSection;
    
    sect->m_dim.clear();
    isOk = m_acdGeomSvc.fillTileData(m_acdId,iVol,sect->m_trans,sect->m_dim,sect->m_center,sect->m_corners);
    if ( ! isOk ) {
      std::cerr << "Failed to fillTileData " << m_acdId.id() << ' ' << iVol << std::endl;
      return StatusCode::FAILURE;
    }
    sect->m_invTrans =  sect->m_trans.inverse();
    if ( sect->m_dim.size() == 5 ) {
      sect->m_trapezoid = true;
      sect->m_xmean = sect->m_dim[0] + sect->m_dim[3];
      sect->m_xmean /= 2.;
      sect->m_xdiff = sect->m_dim[4];      
    }else {
      sect->m_trapezoid = false;
      sect->m_xmean = 0.;
      sect->m_xdiff = 0.;
    }
    m_sections.push_back(sect);
  }
  
  if ( nVol == 2 ) {
    isOk = m_acdGeomSvc.fillTileSharedEdgeData(m_acdId,m_sections[0]->m_dim,m_sections[1]->m_dim,
					       m_sections[0]->m_shared, m_sections[0]->m_shared,
					       m_sections[0]->m_sharedWidth,m_sections[1]->m_sharedWidth);
  }
  
  if ( ! isOk ) return StatusCode::FAILURE;
  return StatusCode::SUCCESS;
}


StatusCode AcdTileDim::latchGaps() {
  
  if ( m_doneGaps ) return StatusCode::SUCCESS;
  if ( m_sections.size() == 0 ) return StatusCode::FAILURE;
  const AcdTileSection* sect = m_sections[0];

  HepPoint3D c1, c2, l1, l2;
  double aX(0.), aY(0.);
  bool isRealGap(false);

  if ( m_acdGeomSvc.getNextTileCorners( m_acdId, -1, c1, c2, isRealGap) != StatusCode::SUCCESS ) {
    return StatusCode::FAILURE;
  }
  if ( isRealGap ) {
    l1 = sect->m_trans * c1;
    l2 = sect->m_trans * c2;
    
    sect->activeDistance(l1,0,aX,aY);
    m_minusXGap[0] = aX < 0 ? -1. * aX : 0;
    sect->activeDistance(l2,0,aX,aY);
    m_minusXGap[1] = aX < 0 ? -1. * aX : 0;
  }    

  if ( m_acdGeomSvc.getNextTileCorners( m_acdId, 1, c1, c2, isRealGap) != StatusCode::SUCCESS ) {
    return StatusCode::FAILURE;
  }
  if ( isRealGap ) { 
    l1 = sect->m_trans * c1;
    l2 = sect->m_trans * c2;
    
    sect->activeDistance(l1,0,aX,aY);
    m_plusXGap[0] = aX < 0 ? -1. * aX : 0;
    sect->activeDistance(l2,0,aX,aY);
    m_plusXGap[1] = aX < 0 ? -1. * aX : 0;
  }

  m_doneGaps = true;
  return StatusCode::SUCCESS;
}


float AcdTileDim::gapAtLocalY(int side, double localY) const {
  float gap = 0.;
  float yFrac = localY / m_sections[0]->m_dim[1];
  yFrac += 0.5;
  float yFrac2 = 1. - yFrac;
  switch ( side ) {
  case -1:
    gap += m_minusXGap[0] * yFrac2;
    gap += m_minusXGap[1] * yFrac;
    break;
  case 1:
    gap += m_plusXGap[0] * yFrac2;
    gap += m_plusXGap[1] * yFrac;
    break;    
  }
  if ( gap > 20. ) {
    std::cout << "Gap at LocalY " << side << ' ' << localY << ' ' << yFrac << ' ' << yFrac2 << std::endl;
    std::cout << m_minusXGap[0] << ' ' << m_minusXGap[1] << ' ' << gap << std::endl;
    std::cout << m_plusXGap[0] << ' ' << m_plusXGap[1] << ' ' << gap << std::endl;
  }
  return gap;
}

float AcdTileDim::gapSize(const HepPoint3D& localPoint, int iVol) const {
  
  if ( iVol == 1 ) {
    return 0.;
  } 
  return gapAtLocalY( (localPoint.x() < 0 ? -1 : 1) , localPoint.y() );
}


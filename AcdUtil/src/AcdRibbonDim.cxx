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

AcdRibbonSegment::AcdRibbonSegment(float halfWidth, 
				   const HepPoint3D& start, const HepPoint3D& end, 
				   const HepTransform3D& trans )
  :m_halfWidth(halfWidth),
   m_start(start),
   m_end(end),
   m_vect(end-start),
   m_trans(trans),
   m_invTrans(trans.inverse()){
}


/// Constructor: takes the ribbon id, the volume id, and the detector service
AcdRibbonDim::AcdRibbonDim(const idents::AcdId& acdId, IAcdGeometrySvc& acdGeomSvc) 
  :m_acdId(acdId),
   m_acdGeomSvc(acdGeomSvc),
   m_halfLength(0.),
   m_topIdx(-1),
   m_plusIdx(-1){
  m_sc = getVals();
}  

/// this function access the detector service to get the geometry information
StatusCode AcdRibbonDim::getVals() {

  StatusCode sc = StatusCode::SUCCESS;

  // First we load up the rays
  bool isOk = m_acdGeomSvc.fillRibbonData(m_acdId,m_segments,m_topIdx,m_plusIdx);
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
  const AcdRibbonSegment* seg = getSegment(isegment);
  if ( seg == 0 ) return;
  local = seg->m_trans * global;
}


double AcdRibbonDim::calcRibbonLength() {
   double len(0);   
   m_segStarts.resize( m_segments.size() + 1 );
   unsigned i(0);
   for ( std::vector<AcdRibbonSegment*>::const_iterator itr = m_segments.begin(); itr != m_segments.end(); itr++, i++ ) {
     float segLen = (*itr)->m_vect.mag();
     m_segStarts[i] = len;
     len += segLen;
   }
   // Total length
   m_segStarts[m_segments.size()] = len;
   return len;  
}

float AcdRibbonDim::getRibbonLength(unsigned iSeg, float segLength) const {
  float retVal(0.);
  if ( (int)iSeg < m_plusIdx ) { 
    // -side or top: Segment frame runs along ribbon, add the segment length to start point for this segment
    retVal = m_segStarts[iSeg] + segLength;
  } else {    
    // +side: Segment frame runs against ribbon, subtract the segment length from start point of next segment
    retVal = m_segStarts[iSeg] + segLength;
  }
  retVal -= m_halfLength;
  return retVal;
}


bool AcdRibbonDim::setEdgeRay(const int& side, HepPoint3D& start, HepVector3D& vect) const{
  HepPoint3D end;
  switch ( side ) {
  case -1:
    start = m_segments[0]->m_start;
    end = m_segments[m_topIdx-1]->m_end;
    break;
  case 0:
    start = m_segments[m_topIdx]->m_start;
    end = m_segments[m_plusIdx-1]->m_end;
    break;
  case 1:
    start = m_segments[m_plusIdx]->m_start;
    end = m_segments.back()->m_end;
    break;
  default:
    return false;
  }
  vect = end - start;
  return true;
}


void AcdRibbonDim::getSegmentsIndices(const HepVector3D& tkDir, bool upward, int& start, int& end, int& dir) const {
  
  switch ( m_acdId.ribbonOrientation() ) {
  case 5:
    dir = tkDir.x() > 0 ? 1 : -1;
    break;
  case 6:
    dir = tkDir.y() > 0 ? 1 : -1;
    break;
  }
  
  switch ( dir ) {
  case -1:
    start = upward ? m_plusIdx - 1 : m_topIdx - 1;
    end = -1;
    break;
  case 1:
    start = upward ? m_topIdx : m_plusIdx;
    end = m_segments.size();
    break;
  }
  

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

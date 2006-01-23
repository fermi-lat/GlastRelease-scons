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

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "CLHEP/Geometry/Transform3D.h"

/// Constructor: takes the ribbon id, the volume id, and the detector service
AcdRibbonDim::AcdRibbonDim(const idents::AcdId& acdId, const idents::VolumeIdentifier& volId, 
			   IGlastDetSvc &detSvc) 
  :m_acdId(acdId),
   m_volId(volId),
   m_detSvc(detSvc) 
{
    m_sc = getVals();

}  

/// this function access the detector service to get the geometry information
StatusCode AcdRibbonDim::getVals() {

  StatusCode sc = StatusCode::SUCCESS;
  static const int ribbonX = 5, ribbonY = 6;
  
  // Need to reconstruct the required volume ids to retrieve the geometry information
  // For now we brute force it.. and construct what we need
  int topSegment;
  int sideFace[2];
  // check orientation to determine which segment number to retrieve for top
  if (m_acdId.ribbonOrientation() == ribbonX) {
    topSegment = 1;
    // ribbons that are along x-axis on the top go down faces 1,3
    sideFace[0] = 1; sideFace[1] = 3;
  } else {
    topSegment = 2;
    // ribbons that are along the y-axis on the top go down faces 2,4
    sideFace[0] = 2; sideFace[1] = 4;
  }
  
  for (int isegment = 0; isegment < 3; isegment++) {
    idents::VolumeIdentifier& segmentVolId = m_segId[isegment];
    segmentVolId.append(m_volId[0]);
    short face;
    if (isegment == 1) {
      face = 0;
    } else {
      face = sideFace[isegment/2];
    }
    segmentVolId.append(face);
    segmentVolId.append(m_volId[2]);
    segmentVolId.append(m_volId[3]);
    segmentVolId.append(m_volId[4]);
    if (isegment != 1) { // side ribbon segment
      segmentVolId.append(0); //put back
    } else { // top ribbon segment
      segmentVolId.append(topSegment);
    }
    
    // Variables for storing ribbon data
    HepPoint3D center;
    double x1=0.0, y1=0.0, z1=0.0;
    double x2=0.0, y2=0.0, z2=0.0;
    double ribbonHalfWidth(0.);
    
    // in this case, we need to extract the dimensions from 3 other top segments
    // to extend an imaginary ribbon across the whole top of the instrument
    if (m_acdId.ribbonOrientation() == ribbonY && isegment == 1) {
      int iseg;
      for(iseg = 1; iseg <= 3; iseg++) {
	idents::VolumeIdentifier volId1;
	volId1.append(m_volId[0]); volId1.append(0); volId1.append(m_volId[2]); 
	volId1.append(m_volId[3]); volId1.append(m_volId[4]);
	if (iseg == 3) {
	  // grab the last segment
	  volId1.append(5);
	} else {
	  volId1.append(iseg);
	}
	
	std::vector<double> dim1;
	sc = getDetectorDimensions(volId1, m_detSvc, dim1, center);
	
	if (sc.isFailure()) {
	  std::cout << "Failed to get dimensions for " << volId1.name() << std::endl;
	   // catch errors
	  sc = StatusCode::SUCCESS;
	  return sc;
	}
	// pick up the beginning from the first segment
	if (iseg == 1) y1 = center.y() - dim1[1]/2.;
	if (iseg == topSegment){
	  // pick up the other 2 dimensions from a ribbon in the middle
	  x1 = center.x(); 
	  z1 = center.z();
	  x2 = center.x();
	  z2 = center.z();
	  ribbonHalfWidth = dim1[0]/2.;
	} else if (iseg == 3) {
	  // Pick up the ending point from the last segment
	  y2 = center.y() + dim1[1]/2.;
	}
      }

    } else if (m_acdId.ribbonOrientation() == ribbonX && isegment == 1) {
      std::vector<double> dim;
      sc = getDetectorDimensions(segmentVolId, m_detSvc, dim, center);
      if (sc.isFailure()) {
	std::cout << "Failed to get dimensions for " <<  segmentVolId.name() << std::endl;
	// catch errors
	sc = StatusCode::SUCCESS;
	continue;
      }
      x1 = center.x() - dim[0]/2.;
      y1 = center.y();
      z1 = center.z();
      x2 = center.x() + dim[0]/2.;
      y2 = center.y();
      z2 = center.z();
      ribbonHalfWidth = dim[1]/2.;

    } else { // side ribbons - which are in one segment
      std::vector<double> dim;
      sc = getDetectorDimensions(segmentVolId, m_detSvc, dim, center);
      if (sc.isFailure()) {
	std::cout << "Failed to get dimensions for " <<  segmentVolId.name() << std::endl;
	// catch errors
	sc = StatusCode::SUCCESS;
	continue;
      }       
      x1 = center.x();
      y1 = center.y();
      z1 = center.z() - dim[2]/2.;
      x2 = center.x();
      y2 = center.y();
      z2 = center.z() + dim[2]/2.;
      // Determine the half width of the ribbon
      if ((face == 1) || (face == 3)) {
	ribbonHalfWidth = dim[0]/2.;
      } else {
	ribbonHalfWidth = dim[0]/2.;
      }
    }
    
    // After all that, we have the beginning and ending points for a line segment
    // that defines a ribbon
    m_start[isegment] = HepPoint3D(x1, y1, z1);
    m_end[isegment] = HepPoint3D(x2, y2, z2);    
    m_halfWidth[isegment] = ribbonHalfWidth;

  }
  
  return sc;
}


StatusCode AcdRibbonDim::getDetectorDimensions( const idents::VolumeIdentifier& volId, IGlastDetSvc &detSvc,
						std::vector<double>& dim, HepPoint3D& center) {

  std::string str;
  StatusCode sc = StatusCode::SUCCESS;
  sc = detSvc.getShapeByID(volId, &str, &dim);
  if ( sc.isFailure() ) {
    //log << MSG::WARNING << "Failed to retrieve Shape by Id" << endreq;
    return sc;
  } 

  HepTransform3D transform;
  sc = detSvc.getTransform3DByID(volId, &transform);
  if (sc.isFailure() ) {
    //log << MSG::WARNING << "Failed to get trasnformation" << endreq;
    return sc;
  } 

  HepPoint3D theCenter(0., 0., 0.);
  center = transform * theCenter;
  return sc;
}

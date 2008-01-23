#include "../AcdRecon/AcdReconFuncs.h"
#include "../AcdRecon/AcdTkrParams.h"

#include "CLHEP/Matrix/Matrix.h"
#include "AcdUtil/AcdTileDim.h"
#include "AcdUtil/AcdRibbonDim.h"
#include "AcdUtil/RayDoca.h"
#include "GlastSvc/Reco/IPropagator.h" 

namespace AcdRecon {

  void pointPoca(const AcdRecon::TrackData& aTrack, const HepPoint3D& point,
		 double& arcLength, double& doca, HepPoint3D& poca) {
    // Poca to a point.  Done global frame

    // dX = vector from the Tile Center to the track start point
    HepVector3D dX = point - aTrack.m_point;

    // arclength to the poca = dot product between track and 
    arcLength = dX * aTrack.m_dir;

    // now calculate the doca    
    doca = sqrt(dX.mag2() - arcLength*arcLength);
    // and the POCA
    poca = aTrack.m_point + arcLength * aTrack.m_dir;   
  }

  void crossesPlane(const AcdRecon::TrackData& aTrack, const HepPoint3D& plane, int iFace, 
		    double& arcLength, HepPoint3D& x_isec) {

    /// This one is used for the intersection w/ the putative ACD volume

    // get stuff from the track
    const HepPoint3D& x0 = aTrack.m_point;
    const HepVector3D& t0 = aTrack.m_dir;

    // Figure out where in the plane of this face the trajectory hits
    if(iFace == 0) {// Top Tile. 
      arcLength = (plane.z()-x0.z())/t0.z();	                
    }
    else if(iFace == 1 || iFace == 3) {// X Side Tile 
      arcLength = (plane.x()-x0.x())/t0.x();
    }
    else if(iFace == 2 || iFace == 4) {// Y Side Tile
      arcLength = (plane.y()-x0.y())/t0.y();
    }
       
    // If arcLength is negative... had to go backwards to hit plane... 
    if( arcLength < 0.) return;
    
    x_isec = aTrack.m_point + arcLength* aTrack.m_dir;
  }


  void crossesPlane(const AcdRecon::TrackData& aTrack,  const HepTransform3D& plane, 
		    double& arcLength, HepPoint3D& x_isec) {

    // This one used the transformation to the plane
 
    // This is done w.r.t. local Z
    const HepPoint3D& x0 = aTrack.m_point;
    const HepVector3D& t0 = aTrack.m_dir;

    // get the point and vector in the plane
    const HepPoint3D& x_p = plane*x0;
    const HepVector3D& t_p = plane*t0;

    // Figure out where the on the trajectory this track hits z = 0 plane
    arcLength = -x_p.z()/t_p.z();

    // If arcLength is negative... had to go backwards to hit plane... 
    if ( arcLength < 0 ) return;    

    // Propagate track along by arcLength to get GLOBAL intersection 
    x_isec = x0 + arcLength * t0;
  }

  void rayDoca(const AcdRecon::TrackData& aTrack, const HepPoint3D&  c1, const HepPoint3D& c2,
	       AcdUtil::RayDoca& rayDoca, double& edgeLen) {

    // Doca between rays (defined by track and two points) done in global frame.
    // Does not account for finite length of ray defined by two points
    const HepPoint3D& x0 = aTrack.m_point;
    const HepVector3D& t0 = aTrack.m_dir;
    Point trackPos(x0.x(), x0.y(), x0.z());
    Vector trackDir(t0.x(), t0.y(), t0.z());
    Ray track(trackPos, trackDir);

    Point pos1(c1.x(), c1.y(), c1.z());      
    Point pos2(c2.x(), c2.y(), c2.z());  
    Vector dir = pos2 - pos1;
    Ray edge(pos1, dir);
    
    rayDoca.recompute(track,edge);
    edgeLen = dir.magnitude();
  }

  void rayDoca_withCorner(const AcdRecon::TrackData& track, const Ray& ray,
			  double& arcLength, double& rayLength, double& dist, Point& x, Vector& v, int& region) {
    // Doca between rays (defined by track and two points)  done in global frame
    // Account for finite length of ray defined by two points.
    // If poca occurs beyond end of that ray, returns poca to end point instead
    const HepPoint3D& x0 = track.m_point;
    const HepVector3D& t0 = track.m_dir;    
    Point trackPos(x0.x(), x0.y(), x0.z());
    Vector trackDir(t0.x(), t0.y(), t0.z());
    Ray trackRay(trackPos, trackDir);
    AcdUtil::RayDoca raydoca(trackRay,ray);
    double alongRay = raydoca.arcLenRay2();
    if ( alongRay < 0 ) {
      HepPoint3D corner(ray.position().x(),ray.position().y(),ray.position().z());
      HepPoint3D poca;
      AcdRecon::pointPoca(track,corner,arcLength,dist,poca);
      x.set(poca.x(),poca.y(),poca.z());
      HepVector3D pocaVect = poca - corner;
      v.set(pocaVect.x(),pocaVect.y(),pocaVect.z());
      rayLength = 0;
      region = -1;
    } else if ( alongRay > ray.getArcLength() ) {
      Point end = ray.position( ray.getArcLength() );
      HepPoint3D corner(end.x(),end.y(),end.z());
      HepPoint3D poca;
      AcdRecon::pointPoca(track,corner,arcLength,dist,poca);
      x.set(poca.x(),poca.y(),poca.z());
      HepVector3D pocaVect = poca - corner;
      v.set(pocaVect.x(),pocaVect.y(),pocaVect.z());
      rayLength = ray.getArcLength();
      region = 1;
    } else {
      arcLength = raydoca.arcLenRay1();
      dist = raydoca.docaRay1Ray2();
      x = raydoca.docaPointRay1();
      v = x - raydoca.docaPointRay2();
      rayLength = alongRay;
      region = 0;
    }       
  }

  void tilePlane(const AcdRecon::TrackData& track, const AcdTileDim& tile,
		 AcdRecon::PocaData& data) {

    // Get all the relevent values for intersections with a tile

    HepPoint3D testHitPoint, testLocal;
    double testActiveX(0.),  testActiveY(0.);
    
    // loop over volumes of tile
    for (int iVol = 0; iVol < tile.nVol(); iVol++ ) {      
      
      double testArcLength(-1.);
      AcdRecon::crossesPlane(track,tile.transform(iVol),testArcLength, testHitPoint);      

      // If arcLength is negative... had to go backwards to hit plane... 
      if ( testArcLength < 0 ) continue;
      
      AcdRecon::tilePlaneActiveDistance(tile,iVol,testHitPoint,testLocal,testActiveX,testActiveY);

      // check to see which active distance is more negative (ie, farther from the center of the tile)
      double testActive2D =  testActiveX < testActiveY ? testActiveX : testActiveY;

      // check to see if the track came closest to this volume
      // if it did, grab the values
      if ( testActive2D > data.m_active2D ) {
	data.m_active2D = testActive2D;
	data.m_activeX = testActiveX;
	data.m_activeY = testActiveY;
	data.m_arcLengthPlane = testArcLength;
	data.m_hitsPlane.set(testHitPoint.x(),testHitPoint.y(),testHitPoint.z());
	data.m_inPlane.set(testLocal.x(),testLocal.y(),testLocal.z());
	data.m_volume = iVol;
      }
    }
  }

  
  void tilePlaneActiveDistance(const AcdTileDim& tile, int iVol, const HepPoint3D& globalPoint,
			       HepPoint3D& localPoint, double& activeX, double& activeY) {

    tile.toLocal(globalPoint,localPoint,iVol);
    tile.activeDistance(localPoint,iVol,activeX,activeY);
  }


  void tileEdgePoca(const AcdRecon::TrackData& aTrack, const AcdTileDim& tile, 
		    double& arcLength, double& dist, Point& x, Vector& v, int& region) {
    // Loop over all edges allowing only the edges which intersect to form
    // the nearest corner participate.
    // For each pair of corners, make a ray and calculate doca from track
    // Note: this could be done at the end - if no edge solution was found however
    //       doing the full monte does not use much more cpu    

    // This is all done in the global frame
    dist = maxDocaValue;
    region = -1;
    
    // loop over volumes of tile
    for (int iVol = 0; iVol < tile.nVol(); iVol++ ) {

      const HepPoint3D* corner = tile.corner(iVol);
      for (int iCorner = 0; iCorner<4; iCorner++) {

	// ignored shared edges
	if ( iCorner == tile.sharedEdge(iVol) ) continue;

	//  WBA: Naively I thought one could limit the edges investigated - not so!
	//	if(iCorner != i_near_corner && (iCorner+1)%4 != i_near_corner) continue;
	int ic2 = iCorner == 3 ? 0 : iCorner + 1;
	const HepPoint3D& c1 = corner[iCorner];
	const HepPoint3D& c2 = corner[ic2];  
	
	// Will need this to determine limit of the tile edge 
	double edge_length(0.);
	
	// Compute DOCA and DOCA location between the track and edge
	AcdUtil::RayDoca raydoca;
	AcdRecon::rayDoca(aTrack,c1,c2,raydoca,edge_length);

	// Make sure that we are going forward
	double s_to_inter = raydoca.arcLenRay1();
	if ( s_to_inter < 0 ) continue;

	// Check if x,y,z along edge falls within limits of tile edge.
	double length_2_intersect = raydoca.arcLenRay2();
	
	if (length_2_intersect > 0 && length_2_intersect < edge_length) {
	  double test_dist = raydoca.docaRay1Ray2();
	  
	  if ( test_dist < dist ) {
	    dist = test_dist;
	    arcLength = raydoca.arcLenRay1();
	    x = raydoca.docaPointRay1();
	    v = x - raydoca.docaPointRay2();
	    region = iVol;
	  }
	}
      }          
    }
  }
  
  void tileEdgeCornerPoca(const AcdRecon::TrackData& aTrack, const AcdTileDim& tile, 
			  double& arcLength, double& dist, Point& x, Vector& v, int& region) {

    // Get four corners associated with the tile.
    // Assuming we can avoid, calculation with all 8 corners, 4 should be enough
    // The corners are returned in order (-,-), (-,+), (+,+), (+,-)
    // where third dimension is the one we ignore, since it is associated with
    // tile thickness.

    // This is all done in the global frame
    region = -1;
    dist = maxDocaValue;
    
    // loop over volumes of tile
    for (int iVol = 0; iVol < tile.nVol(); iVol++ ) {

      const HepPoint3D* corner = tile.corner(iVol);

      // First find the nearest corner and the distance to it. 
      double arclen(-1.);    
      double test_dist(2000);
      for (int iCorner = 0; iCorner<4; iCorner++) {      
	HepPoint3D x_isec;
	test_dist = 2000.;
	AcdRecon::pointPoca(aTrack,corner[iCorner],arclen,test_dist,x_isec);
	if ( test_dist < dist ) {
	  dist = test_dist;
	  region = iVol;
	  arcLength = arclen;
	  x.set(x_isec.x(),x_isec.y(),x_isec.z());
	  HepVector3D pocaVect = x_isec - corner[iCorner];
	  v.set(pocaVect.x(),pocaVect.y(),pocaVect.z());
	}
      }    
      Point pEdge;
      Vector vEdge;
      test_dist = 2000.;
      int test_region(-1);
      AcdRecon::tileEdgePoca(aTrack,tile,arclen,test_dist,pEdge,vEdge,test_region);
      if ( test_dist < dist ) {
	dist = test_dist;
	region = test_region;
	arcLength = arclen;
	x = pEdge;
	v = vEdge;
      }
    }
    // Make this an Active Distance calculation (we know that we are outside the tile)
    dist *= -1.;
  }

  void ribbonPlane(const AcdRecon::TrackData& aTrack, const AcdRibbonDim& ribbon,
		   AcdRecon::PocaData& data){
    
    // FIXME, use the transformations instead

    // get the ribbon id
    // const idents::AcdId& acdId = ribbon.acdId();  
    
    // this is the value we want to beat
    double best_dist = -maxDocaValue;
    
    HepPoint3D testHitPoint, testLocal;

    // loop over segments
    for ( int segment = 0; segment < 3; segment++ ) {

      double testArcLength(-1.);
      const HepTransform3D& toFacePlane = ribbon.transform(segment);
      AcdRecon::crossesPlane(aTrack, toFacePlane, testArcLength, testHitPoint);      
	
      if (testArcLength < 0) continue;

      testLocal = toFacePlane * testHitPoint;
    
      double test_dist = fabs(testLocal.x());
      
      // Make this an Active Distance calculation 
      test_dist = ribbon.halfWidth() - test_dist;

      // test against the best value
      if ( test_dist > best_dist ) {
	data.m_active2D = test_dist;
	best_dist = test_dist;
	data.m_arcLengthPlane = testArcLength;
	data.m_hitsPlane.set(testHitPoint.x(),testHitPoint.y(),testHitPoint.z());
	data.m_volume = segment;
      }
    }
  }


  void ribbonPoca(const AcdRecon::TrackData& track, const AcdRibbonDim& ribbon,
		  double& arcLength, double& ribbonLength, double& dist, Point& x, Vector& v, int& region) {

    // 
   
    // First, does this ribbon extend along x or y axix
    static const int ribbonX = 5;
    const idents::AcdId& acdId = ribbon.acdId();  
    int dir = 0;
    if ( acdId.ribbonOrientation() == ribbonX ) {
      // extends along x.  Want to know if it is going towards +-Y sides  
      dir = track.m_dir.x() > 0 ? 1 : -1;
    } else {
      // extends along y.  Want to know if it is going towards +-X sides  
      dir = track.m_dir.y() > 0 ? 1 : -1;
    }

    dist = -2000.;
    double arcLengthTest(0.);
    double distTest(0.);
    Point x_test;
    Vector v_test;
    int regionTest(0);
    double dist_last(-500000.);
    std::vector<const Ray*> raysInOrder;    
    
    int iRay(0);
    
    // Now, look over the relevent rays.
    if ( dir == 1 ) {
      // Going to + side
      for ( iRay = 0; iRay < (int)ribbon.plusSideRays().size(); iRay++ ) {
	const Ray& aRay = ribbon.plusSideRays()[iRay];
	raysInOrder.push_back(&aRay);
      }
      if ( track.m_upward ) {
	for ( iRay = ribbon.topRays().size() -1; iRay >= 0; iRay-- ) {
	  const Ray& aRay = ribbon.topRays()[iRay];
	  raysInOrder.push_back(&aRay);
	}
      }
    } else if ( dir == -1 ) {
      // Going to - side
       for ( iRay = 0; iRay < (int)ribbon.minusSideRays().size(); iRay++ ) {
	const Ray& aRay = ribbon.minusSideRays()[iRay];
	raysInOrder.push_back(&aRay);
       }
       if ( track.m_upward ) {
	 for ( iRay =0; iRay < (int)ribbon.topRays().size(); iRay++ ) {
	   const Ray& aRay = ribbon.topRays()[iRay];
	   raysInOrder.push_back(&aRay);
	 }
       }
    } else {
      return;
    }

    double ribbonLen(0);
    double rayLen(0);

    for ( iRay = 0; iRay < (int)raysInOrder.size(); iRay++ ) {
      const Ray& aRay = *(raysInOrder[iRay]);
      AcdRecon::rayDoca_withCorner(track,aRay,arcLengthTest,rayLen,distTest,x_test,v_test,regionTest);      

      // flip the ray length to top rays and +x,+y going.
      if ( dir == 1 && iRay >= (int)ribbon.plusSideRays().size() ) {
	rayLen = aRay.getArcLength() - rayLen;
      }

      // Make this an Active Distance calculation 
      distTest = ribbon.halfWidth() - distTest;

      if ( distTest < dist_last ) {
	// going wrong direction. stop
	//return;
	;
      }
      dist_last = distTest;
      if ( distTest >  dist ) {
	dist = distTest;
	arcLength = arcLengthTest;
	ribbonLength = ribbon.halfLength() - (ribbonLen + rayLen  );
	ribbonLength *= dir;
	x = x_test;
	v = v_test;
	region = regionTest;
      }
      ribbonLen += aRay.getArcLength();      
    }
  }

  // initiliaze the kalman propagator
  void startPropagator(IPropagator& prop, const Event::TkrTrack& aTrack, const AcdRecon::TrackData& trackData,
		       const double& maxArcLength ) {

    // Get the start point & direction of the track & the params & energy also
    const unsigned int hitIndex = trackData.m_upward ? 0 : aTrack.getNumHits() - 1;
    const Event::TkrTrackHit* theHit = aTrack[hitIndex];
    const Point initialPosition = theHit->getPoint(Event::TkrTrackHit::SMOOTHED);
    const Event::TkrTrackParams& trackPars = theHit->getTrackParams(Event::TkrTrackHit::SMOOTHED); 
    
    // setup the propagator
    prop.setStepStart(trackPars,initialPosition.z(),trackData.m_upward); 
    prop.step(maxArcLength);  
  }


  // run the propagtor out to a specified arclength
  void propagateToArcLength(IPropagator& prop,
			    const AcdRecon::TrackData& track, const double& arcLength,
			    Event::TkrTrackParams& next_params,
			    AcdTkrParams& paramsAtArcLength) {

    Point x_step = prop.getPosition( arcLength );
    next_params = prop.getTrackParams(arcLength,track.m_energy,true);
    paramsAtArcLength.set(next_params,x_step.z(),track.m_upward);
  }


  void projectErrorToPlane(const AcdTkrParams& paramsAtArcLength, const CLHEP::HepMatrix& localFrameVectors,
			  CLHEP::HepSymMatrix& covAtPlane) {
    //  U = A V A^T
    //

    for ( int i(1); i < 3; i++ ) {
      for ( int j = i; j < 3; j++) {
	double currentSum = 0.;
	for ( unsigned int k(1);  k < 4; k++ ) {
	  for ( unsigned int l(1); l < 4; l++ ) {
	    currentSum += localFrameVectors(i,k) * paramsAtArcLength(k,l) * localFrameVectors(j,l);
	  }
	}	
	covAtPlane(i,j) = currentSum;
      }
    }
  }

  
  // Project error ellipse onto a line
  void projectErrorToPocaVector(const AcdTkrParams& paramsAtArcLength, const Vector& pocaVector, 
				double& pocaError) {
    //  U = A V A^T
    //
    pocaError = 0.;
    for ( unsigned int i(0);  i < 3; i++ ) {
      for ( unsigned int j(0);  j < 3; j++ ) {
	pocaError += pocaVector[i] * paramsAtArcLength(i+1,j+1) * pocaVector[j];	  
      }
    }
  }

  bool exitsLat(const AcdRecon::TrackData& trackData,
		const AcdRecon::AcdVolume& acdVol,
		AcdRecon::ExitData& data) {

    // grab the track data
    const HepPoint3D& initialPosition = trackData.m_point;
    const HepVector3D& initialDirection = trackData.m_dir;    

    // sanity check
    if ( (trackData.m_upward && initialDirection.z() < 0) || (!trackData.m_upward && initialDirection.z() > 0) ) {
      return false;
    }

    // hits -x or +x side ?
    const double normToXIntersection =  initialDirection.x() < 0 ?  
      -1.*acdVol.m_sides - initialPosition.x() :    // hits -x side
      1.*acdVol.m_sides - initialPosition.x();      // hits +x side  
    const double slopeToXIntersection = fabs(initialDirection.x()) > 1e-9 ? 
      1. / initialDirection.x() : (normToXIntersection > 0. ? 1e9 : -1e9);
    const double sToXIntersection = normToXIntersection * slopeToXIntersection;
    
    // hits -y or +y side ?
    const double normToYIntersection = initialDirection.y() < 0 ?  
      -1.*acdVol.m_sides - initialPosition.y() :    // hits -y side
      1.*acdVol.m_sides - initialPosition.y();      // hits +y side
    const double slopeToYIntersection = fabs(initialDirection.y()) > 1e-9 ? 
      1. / initialDirection.y() : (normToYIntersection > 0. ? 1e9 : -1e9); 
    const double sToYIntersection = normToYIntersection * slopeToYIntersection;
    
    // hits top or bottom
    const double normToZIntersection =  trackData.m_upward ? 
      acdVol.m_top - initialPosition.z() :
      acdVol.m_bottom - initialPosition.z();
    const double slopeToZIntersection = 1. / initialDirection.z();
    const double sToZIntersection = normToZIntersection * slopeToZIntersection;    

    // pick the closest plane
    if ( sToXIntersection < sToYIntersection ) {
      if ( sToXIntersection < sToZIntersection ) {
	// hits X side
	data.m_arcLength = sToXIntersection;
	data.m_face = initialDirection.x() > 0 ? 1 : 3;
      } else {
	// hits Z side
	data.m_arcLength = sToZIntersection;
	data.m_face = trackData.m_upward ? 0 : 5;
      }     
    } else {
      if ( sToYIntersection < sToZIntersection ) {
	// hits Y side
	data.m_arcLength = sToYIntersection;
	data.m_face = initialDirection.y() > 0 ? 2 : 4;
      } else {
	// hits Z side
	data.m_arcLength = sToZIntersection;
	data.m_face =  trackData.m_upward ? 0 : 5;
      }     
    }
    
    // protect against negative arcLengths
    if ( data.m_arcLength < 0. ) {
      return false;
    }
    
    // extrapolate to the i-sect
    HepPoint3D iSect = initialPosition + data.m_arcLength*initialDirection;
    data.m_x.set(iSect.x(),iSect.y(),iSect.z());
    
    // flip the sign of the arclength for downgoing side
    data.m_arcLength *= trackData.m_upward ? 1. : -1.;
    
    return true;
  }

  bool entersLat(const AcdRecon::TrackData& trackData, 
		 const AcdRecon::AcdVolume& acdVol,
		 AcdRecon::ExitData& data) {
    
    // grab the track data
    const HepPoint3D& initialPosition = trackData.m_point;
    const HepVector3D& initialDirection = trackData.m_dir;    
           
    bool enters(false);

    // where does the track start relative to +-X sides
    // and how long before it hits one of the sides
    // this evals to -1 if between sides
    double sToXIntersection(-1.);
    if ( fabs(initialPosition.x()) > acdVol.m_sides ) {
      double normToXIntersection = initialPosition.x() < 0 ? 
	-1.*acdVol.m_sides - initialPosition.x() :    // hits -x side first
	1.*acdVol.m_sides - initialPosition.x();      // hits +x side frist
      const double slopeToXIntersection = fabs(initialDirection.x()) > 1e-9 ? 
	1. / initialDirection.x() : (normToXIntersection > 0. ? 1e9 : -1e9);
      sToXIntersection = normToXIntersection * slopeToXIntersection;
      // propagate to that point, make sure that other two values inside LAT also
      if ( sToXIntersection > 0 ) {
	HepPoint3D xPlaneInter = initialPosition;  xPlaneInter += sToXIntersection* initialDirection;
	if (  fabs(xPlaneInter.y()) < acdVol.m_sides  &&
	      xPlaneInter.z() > acdVol.m_bottom &&
	      xPlaneInter.z() < acdVol.m_top ) {
	  data.m_arcLength = sToXIntersection;
	  data.m_x.set(xPlaneInter.x(),xPlaneInter.y(),xPlaneInter.z());
	  data.m_face = initialPosition.x() < 0 ? 1 : 3;
	  enters = true;
	}
      }
    }  

    // where does the track start relative to +-Y sides
    // this evals to -1 if between sides
    double sToYIntersection(-1.);
    if ( fabs(initialPosition.y()) > acdVol.m_sides ) {
      double normToYIntersection = initialPosition.y() < 0 ? 
	-1.*acdVol.m_sides - initialPosition.y() :    // hits -y side first
	1.*acdVol.m_sides - initialPosition.y();      // hits +y side frist
      const double slopeToYIntersection = fabs(initialDirection.y()) > 1e-9 ? 
	1. / initialDirection.y() : (normToYIntersection > 0. ? 1e9 : -1e9);
      sToYIntersection = normToYIntersection * slopeToYIntersection;
      if ( sToYIntersection > 0 && 
	   ( sToYIntersection < data.m_arcLength || data.m_arcLength < 0 ) ) {    
	// propagate to that point, make sure that other two values inside LAT also
	HepPoint3D yPlaneInter = initialPosition;  yPlaneInter += sToYIntersection* initialDirection;
	if (  fabs(yPlaneInter.x()) < acdVol.m_sides  &&
	      yPlaneInter.z() > acdVol.m_bottom &&
	      yPlaneInter.z() < acdVol.m_top ) {
	  data.m_arcLength = sToYIntersection;
	  data.m_x.set(yPlaneInter.x(),yPlaneInter.y(),yPlaneInter.z());
	  data.m_face = initialPosition.y() < 0 ? 2 : 4;
	  enters = true;
	}
      }
    }  

    // where does the track start relative to +-Z sides
    // this evals to -1 if between sides
    double sToZIntersection(-1.);
    if ( initialPosition.z() < acdVol.m_bottom ||
	 initialPosition.z() > acdVol.m_top ) {
      double normToZIntersection = initialPosition.z() < 0 ? 
	acdVol.m_bottom - initialPosition.z() :    // hits -z side first
	acdVol.m_top - initialPosition.z();      // hits +z side frist
      const double slopeToZIntersection = fabs(initialDirection.z()) > 1e-9 ? 
	1. / initialDirection.z() : (normToZIntersection > 0. ? 1e9 : -1e9);
      sToZIntersection = normToZIntersection * slopeToZIntersection;
      if ( sToZIntersection > 0 && 
	   ( sToZIntersection < data.m_arcLength || data.m_arcLength < 0 ) ) {    
	// propagate to that point, make sure that other two values inside LAT also
	HepPoint3D zPlaneInter = initialPosition;  zPlaneInter += sToZIntersection* initialDirection;
	if (  fabs(zPlaneInter.x()) < acdVol.m_sides  &&
	      fabs(zPlaneInter.y()) < acdVol.m_sides ) {
	  data.m_arcLength = sToZIntersection;
	  data.m_x.set(zPlaneInter.x(),zPlaneInter.y(),zPlaneInter.z());
	  data.m_face = initialPosition.z() > acdVol.m_top ? 0 : 5;	
	  enters = true;
	}
      }
    }  
    return enters;
  }
 

  void entersCal(const AcdRecon::TrackData& aTrack, const double& calZDist, 
		 Point& entryPoint, Vector& entryVector, int& region) {

    const HepPoint3D& x0 = aTrack.m_point;
    const HepVector3D& t0 = aTrack.m_dir;

    // set a default
    region = -1;

    // Figure out where in the plane of this face the trajectory hits
    double arcLength = (calZDist - x0.z())/t0.z();	                

    // If arcLength is negative... had to go backwards to hit plane... 
    if( arcLength < 0.) return;

    HepPoint3D x_isec = aTrack.m_point + arcLength* aTrack.m_dir;
    entryPoint.set( x_isec.x(), x_isec.y(), x_isec.z() );
    entryVector.set( t0.x(), t0.y(), t0.z() );

    region = 0;
    if ( entryPoint.x() < -800. ) region += 1;
    if ( entryPoint.x() > 800. ) region += 2;    
    if ( entryPoint.y() < -800. ) region += 4;
    if ( entryPoint.y() > 800. ) region += 8;

    return;
  }

  void splashVariables(const AcdTileDim& tile, const Point& entryPoint, const Vector& entryVector,
		       double& solidAngle, double& meanAngle, double& pathLength) {

    // This code is from Phillipe Bruel, 
    // modified by EAC to match the format of the data in AcdTileDim and to use CLHEP
    //
    // The idea is to split a tile up into many little elements and calculate three quantities for each
    // element and average or sum these quantities over the whole tile.
    //
    //  solidAngle -> solidAngle as seen from the impact point int the CAL
    //  meanAngle -> the angle between the track direction and the vector from impact to tile, weighted by solid angle
    //  pathLength -> the pathLength through the tile, weighted by solid angle
    //
    // 
    //  Stuff used in the calculation
    //  
    // the center of the tile : pc[3]  -> called "center" here
    // 
    // pcal[3] = the impact point of the beam particle onto the top of the calorimeter
    //               called "impact" here
    //
    // velec[3] = the direction of the beam particle (the norm of this vector = 1)
    //               called "direction" here 
    //
    // three vectors : v0[4], v1[4], v2[4]
    // - v0[0,1,2] define the axis, v0[3] = the dimension along this axis (and the same for v1 and v2)
    // - v2 is always the axis corresponding to the smallest dimension (i.e the width of the tile = 10mm)
    // - v2 points towards the outside world : a backsplash particle, coming 
    // from the impact point of the beam particle onto the top of the 
    // calorimeter, would pass through the tile with a momentum vector p that would gives p.v2>0
    // - (v0,v1,v2) is a direct system
    //
    // this logic is just handled by the switch statement below, the equivalent information is in
    // 
    // "stepVector1", "stepVector2" and "normalVector" 
    // which are the vector to step from one element to the next and the unit normal to the tile surface
    //
    //

    // reset are the variables that we want to determine
    solidAngle = 0.;
    meanAngle = 0.;
    pathLength = 0.;
    double ANGLERMS = 0;   // this is getting dropped might want to save it later
    
    // definition of the surface elements steps
    static const int NSTEP = 50;
    static const int nxstep = NSTEP;
    static const int nystep = NSTEP;
    static const int nzstep = NSTEP;

    // the incoming track
    const HepPoint3D impact(entryPoint.x(),entryPoint.y(),entryPoint.z());
    const HepVector3D direction(entryVector.x(),entryVector.y(),entryVector.z());

    // get the tile center and corner
    const HepPoint3D& center = tile.tileCenter();
    const HepPoint3D& lowCorner = (tile.corner(0))[0];

    // define the step vectors
    double halfX = (center.x() - lowCorner.x()) / (double)nxstep;
    double halfY = (center.y() - lowCorner.y()) / (double)nystep;
    double halfZ = (center.z() - lowCorner.z()) / (double)nzstep;    

    HepVector3D halfStepVector;    
    HepVector3D stepVector1;
    HepVector3D stepVector2;
    HepVector3D normalVector;
    int nstep1;
    int nstep2;

    double surfaceElementArea(0.);
    double width(0.);
    switch ( tile.acdId().face() ) {
    case 0:
      // top.  1 -> X, 2 -> Y, thickness -> Z
      halfStepVector = HepVector3D(halfX,halfY,0.);
      stepVector1 = HepVector3D(halfX*2.,0.,0.);
      stepVector2 = HepVector3D(0.,halfY*2.,0.);
      normalVector = HepVector3D(0.,0.,1.);
      surfaceElementArea = halfX * halfY * 4.;
      width = tile.dim()[2];
      nstep1 = nxstep;
      nstep2 = nystep;
      break;
    case 1:
    case 3:
      // +-X.  1 -> Y, 2 -> Z, thickness -> X
      halfStepVector = HepVector3D(0.,halfY,halfZ);
      stepVector1 = HepVector3D(0.,halfY*2.,0.);
      stepVector2 = HepVector3D(0.,0.,halfZ*2.);
      normalVector = HepVector3D(1.,0.,0.);
      surfaceElementArea = halfY * halfZ * 4.;      
      width = tile.dim()[0];
      nstep1 = nystep;
      nstep2 = nzstep;
      break;
    case 2:
    case 4:
      // +-Y.  1 -> X, 2 -> Z, thickness -> Y
      halfStepVector = HepVector3D(halfX,0.,halfZ);
      stepVector1 = HepVector3D(halfX*2.,0.,0.);
      stepVector2 = HepVector3D(0.,0.,halfZ*2.);
      normalVector = HepVector3D(0.,1.,0.);
      surfaceElementArea = halfX * halfZ * 4.;
      width = tile.dim()[1];
      nstep1 = nxstep;
      nstep2 = nzstep;
      break;
    }
    
    switch ( tile.acdId().face() ) {
    case 1:
    case 2:
      // -X or -Y side.  Invert the normal vector.
      normalVector *= -1.;
      break;
    }

    // loop on the surface elements, start a half step inside the 
    HepPoint3D currentSurfacePoint = lowCorner + halfStepVector;
    int iPoint(0);
    for(int i=1;i<nstep1;++i) {
      for(int j=1;j<nstep2;++j) {

	iPoint++;

	// determine the vector (impact point, center of the surface element)
	HepVector3D toElement = currentSurfacePoint - impact;

	double distance_to_tile_squared = toElement.mag2();
	HepVector3D unitToElement = toElement.unit();

	double tiltAngle = unitToElement.dot(normalVector);

	// calculate the solid angle of the surface element seen from the impact point
	double solang = tiltAngle * surfaceElementArea / distance_to_tile_squared;
	
	// add it to the total solidangle
	solidAngle += solang;

	// calculate the angle between the beam particle direction and the 
	// vector (impact point, center of the surface element)
	double flightAngle = acos( unitToElement.dot(direction) );

	// use it to calculate the mean and the rms of the angle with the weight solang
	meanAngle += flightAngle*solang;
	ANGLERMS += flightAngle*flightAngle*solang;
	
	// calculate the corrected width.
	double widthcorrection = unitToElement.dot(normalVector);	
	if(widthcorrection>0)
	  pathLength += (width/widthcorrection)*solang;

	// debugging printout
	//std::cout << currentSurfacePoint << ' ' << impact << ' ' 
	//	  << toElement << ' ' << unitToElement << ' ' 
	//	  << solang << ' ' << flightAngle << ' ' << widthcorrection << std::endl;	
	
	// step in local Y
	currentSurfacePoint += stepVector2;
      }
      // reset localY, and step in local X
      currentSurfacePoint -= ((float)(nstep2-1))*stepVector2;
      currentSurfacePoint += stepVector1;
    }

    float meanSolid = solidAngle / (float)(iPoint);

    // divide by the total weight = the total solid angle
    if(solidAngle>0) {
      meanAngle /= solidAngle;
      ANGLERMS /= solidAngle;
      ANGLERMS = sqrt(ANGLERMS-meanAngle*meanAngle);
      pathLength /= solidAngle;
    }    
    return;
  }

}

#include "../AcdRecon/AcdReconFuncs.h"

#include "CLHEP/Matrix/Matrix.h"

namespace AcdRecon {

  void pointPoca(const AcdRecon::TrackData& aTrack, const HepPoint3D& point,
		 double& arcLength, double& doca, HepPoint3D& poca) {

    // dX = vector from the Tile Center to the track start point
    HepVector3D dX = point - aTrack.m_point;

    // arclength to the poca = dot product between track and 
    arcLength = dX * aTrack.m_dir;

    // If arcLength is negative... had to go backwards to hit plane... 
    // if (  arcLength < 0. ) return;

    // now calculate the doca    
    doca = sqrt(dX.mag2() - arcLength*arcLength);
    // and the POCA
    poca = aTrack.m_point + arcLength * aTrack.m_dir;   
  }

  void crossesPlane(const AcdRecon::TrackData& aTrack, const HepPoint3D& plane, int iFace, 
		    double& arcLength, double& localX, double& localY, HepPoint3D& x_isec) {
 
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
    HepVector3D local_x0 = x_isec - plane;
    
    if(iFace == 0) {// Top Tile
      localX = local_x0.x();
      localY = local_x0.y();
      // Choose which is furthest away from edge (edge @ 0.)
    } else if(iFace == 1 || iFace == 3) {// X Side Tile
      localX = local_x0.y();
      localY = local_x0.z();
    } else if(iFace == 2 || iFace == 4) {// Y Side Tile
      localX = local_x0.x();
      localY = local_x0.z();
    }
  }
  
  
  void rayDoca(const AcdRecon::TrackData& aTrack, const HepPoint3D&  c1, const HepPoint3D& c2,
	       RayDoca& rayDoca, double& edgeLen) {
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
    const HepPoint3D& x0 = track.m_point;
    const HepVector3D& t0 = track.m_dir;    
    Point trackPos(x0.x(), x0.y(), x0.z());
    Vector trackDir(t0.x(), t0.y(), t0.z());
    Ray trackRay(trackPos, trackDir);
    RayDoca raydoca(trackRay,ray);
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
		 double& arcLength, double& localX, double& localY, 
		 double& activeX, double& activeY, double& active2D, HepPoint3D& hitPoint) {

    // loop over volumes of tile
    for (int iVol = 0; iVol < tile.nVol(); iVol++ ) {      
      
      double testActiveX(-maxDocaValue);
      double testActiveY(-maxDocaValue);
      double testLocalX(maxDocaValue);
      double testLocalY(maxDocaValue);
      double testArcLength(-1.);
      HepPoint3D testHitPoint;
      
      std::vector<double> dim = tile.dim(iVol);
      // Beware: these dimensions are in some sort of local system and for
      // iFace = 1 || 3  x<->y 
      double dX = dim[0];
      double dY = dim[1];
      double dZ = dim[2];
      
      const HepPoint3D& center = tile.tileCenter(iVol);

      // fix this
      int face = tile.face(iVol);

      crossesPlane(track,center,face,testArcLength,testLocalX,testLocalY,testHitPoint);

      // int region(0);
 
      // If arcLength is negative... had to go backwards to hit plane... 
      if ( testArcLength < 0 ) continue;
      if(face == 0) {// Top Tile
	testActiveX = dX/2. - fabs(testLocalX);
	testActiveY = dY/2. - fabs(testLocalY);
	// check to see if this is one of the top-side tiles
	if ( tile.nVol() == 2 ) {
	  if ( tile.sharedEdge(0) == 1 && testLocalY > 0 ) {
	    testActiveY += fabs(tile.sharedWidth(0));
	  } else if (  tile.sharedEdge(0) == 3 && testLocalY < 0 ) {
	    testActiveY += fabs(tile.sharedWidth(0));
	  }	  
	}
      } else if(face == 1 || face == 3) {// X Side Tile
	testActiveY = dZ/2  - fabs(testLocalY);
	testActiveX = dY/2. - fabs(testLocalX);
      } else if(face == 2 || face == 4) {// Y Side Tile
	testActiveY = dZ/2. - fabs(testLocalY);
	testActiveX = dX/2. - fabs(testLocalX);
	// check to see if this is one of the extra pieces of the side tiles
	if ( tile.nVol() == 2 ) {
	  if ( tile.sharedEdge(1) == 1 && testLocalY > 0 ) {
	    // is a shared piece.  but this is a short side, so take the distance to the
	    // other side of this volume
	    testActiveY =  dZ - testActiveY;
	    // testActiveY += fabs(tile.sharedWidth(1));
	  }	  
	}
      }

     

      // check to see which active distance is more negative (ie, farther from the center of the tile)
      double testActive2D =  testActiveX < testActiveY ? testActiveX : testActiveY;

      // check to see if the track came closest to this volume
      // if it did, grab the values
      if ( testActive2D > active2D ) {
	active2D = testActive2D;
	activeX = testActiveX;
	activeY = testActiveY;
	localX = testLocalX;
	localY = testLocalY;
	arcLength = testArcLength;
	hitPoint = testHitPoint;
      }
    }
  }

  void tileEdgePoca(const AcdRecon::TrackData& aTrack, const AcdTileDim& tile, 
		    double& arcLength, double& dist, Point& x, Vector& v, int& region) {
    // Loop over all edges allowing only the edges which intersect to form
    // the nearest corner participate.
    // For each pair of corners, make a ray and calculate doca from track
    // Note: this could be done at the end - if no edge solution was found however
    //       doing the full monte does not use much more cpu

    
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
	RayDoca raydoca;
	AcdRecon::rayDoca(aTrack,c1,c2,raydoca,edge_length);

	// Make sure that we are going forward
	double s_to_inter = raydoca.arcLenRay1();
	if ( s_to_inter < 0 ) continue;

	// Check if x,y,z along edge falls within limits of tile edge.
	double length_2_intersect = raydoca.arcLenRay2();
	
	if (length_2_intersect > 0 && length_2_intersect < edge_length) {
	  double test_dist = raydoca.docaRay1Ray2();

	  // check to see if this is a shared edge of a curved tile
	  if ( test_dist < dist ) {
	    dist = test_dist;
	    arcLength = raydoca.arcLenRay1();
	    region = iCorner + 10*iVol;
	    x = raydoca.docaPointRay1();
	    v = x - raydoca.docaPointRay2();
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
	  region = iCorner+4 +(10*iVol);
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
    // Make this an Active Distance calculation 
    dist *= -1.;
  }

  void ribbonPlane(const AcdRecon::TrackData& aTrack, const AcdRibbonDim& ribbon,
		  double& arcLength, double& dist, HepPoint3D& x) {
    
    // get the ribbon id
    const idents::AcdId& acdId = ribbon.acdId();  
    
    // this is the value we want to beat
    double best_dist = -maxDocaValue;

    // loop over segments
    for ( int segment = 0; segment < 3; segment++ ) {

      HepPoint3D sideCenter;
      const HepPoint3D origin;
      ribbon.toLocal(origin,segment,sideCenter);

      // check orientation to determine which segment corresponds to which face
      int face(-1);      
      switch ( segment ) {
      case 1: face = 0; break;
      case 0: face = acdId.ribbonOrientation() == ribbonX ? 1 : 2; break;
      case 2: face = acdId.ribbonOrientation() == ribbonX ? 3 : 4; break;
      }
      
      // where does this track cross the plane of the tile 
      HepPoint3D x_isec;
      HepPoint3D ribbonStartPos;
      HepVector3D ribbonVec;
      double test_arc(-1.);
      double localX(0.); double localY(0.);
      AcdRecon::crossesPlane(aTrack,sideCenter,face,test_arc,localX,localY,x_isec);
	
      if (test_arc < 0) continue;
    
      bool isOk = ribbon.setEdgeRay(segment,ribbonStartPos,ribbonVec);
      if ( !isOk ) continue;
      
      // Form vector between the beginning of the ribbon and the point where the
      // track intersects the plane of the ribbon
      HepVector3D delta = x_isec - ribbonStartPos;
      // Form a vector for the ribbon
      double prod = delta * ribbonVec.unit();
      // check that the projection of the point to the ribbon occurs within the
      // length of the ribbon segment
      if ((prod < 0) || (prod > ribbonVec.mag())) continue;
    
      double test_dist = sqrt(delta.mag2() - prod*prod);

      // Make this an Active Distance calculation 
      test_dist = ribbon.halfWidth() - test_dist;

      // test against the best value
      if ( test_dist > best_dist ) {
	dist = test_dist;
	best_dist = test_dist;
	arcLength = test_arc;
	x = x_isec;
      }
    }
  }


  void ribbonPoca(const AcdRecon::TrackData& track, const AcdRibbonDim& ribbon,
		  double& arcLength, double& ribbonLength, double& dist, Point& x, Vector& v, int& region) {
   
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
      for ( iRay = 0; iRay < ribbon.plusSideRays().size(); iRay++ ) {
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
       for ( iRay = 0; iRay < ribbon.minusSideRays().size(); iRay++ ) {
	const Ray& aRay = ribbon.minusSideRays()[iRay];
	raysInOrder.push_back(&aRay);
       }
       if ( track.m_upward ) {
	 for ( iRay =0; iRay < ribbon.topRays().size(); iRay++ ) {
	   const Ray& aRay = ribbon.topRays()[iRay];
	   raysInOrder.push_back(&aRay);
	 }
       }
    } else {
      return;
    }

    double ribbonLen(0);
    double rayLen(0);

    for ( iRay = 0; iRay < raysInOrder.size(); iRay++ ) {
      const Ray& aRay = *(raysInOrder[iRay]);
      AcdRecon::rayDoca_withCorner(track,aRay,arcLengthTest,rayLen,distTest,x_test,v_test,regionTest);      

      // flip the ray length to top rays and +x,+y going.
      if ( dir == 1 && iRay >= ribbon.plusSideRays().size() ) {
	rayLen = aRay.getArcLength() - rayLen;
      }

      // Make this an Active Distance calculation 
      distTest = ribbon.halfWidth() - distTest;

      if ( distTest < dist_last ) {
	// going wrong direction. stop
	//return;
	continue;
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

  void projectToPlane(const AcdRecon::TrackData& track, const Event::TkrTrackParams& params, 
		      int face, const HepPoint3D& plane,
		      AcdRecon::PocaData& data) {
    
    double delta(0.);
    HepMatrix covAtPlane(2,2);
    switch ( face ) {
    case 0:
      // top: x,y,z are same for local and global
      delta = plane.z() - track.m_point.z();
      data.m_cosTheta = track.m_dir.z();
      AcdRecon::errorAtZPlane(delta,params,covAtPlane);
    break;
    case 1:
    case 3:
      // +-x sides: y -> local x, z -> local y
      delta = plane.x() - track.m_point.x();
      data.m_cosTheta = fabs(track.m_dir.x());
      AcdRecon::errorAtXPlane(delta,params,covAtPlane);
      break;
    case 2:
    case 4:
      // +-y sides: x -> local x, z -> local y
      delta = plane.y() - track.m_point.y();
      data.m_cosTheta = fabs(track.m_dir.y());
      AcdRecon::errorAtYPlane(delta,params,covAtPlane);
      break;
    } 
    data.m_localCovXX = covAtPlane[0][0];
    data.m_localCovYY = covAtPlane[1][1];
    data.m_localCovXY = covAtPlane[0][1];
  }

  
  void errorAtXPlane(const double delta, const Event::TkrTrackParams& pars, HepMatrix& covAtPlane) {
    
    // get the tk params.  
    // want the x and y and the normal to the plane
    const double m_x = pars.getxSlope();
    const double inv_m_x = fabs(m_x) > 1e-9 ? 1./ m_x : 1e9;
    /* not used */ //const double m_y = pars.getySlope();
    
    // ok input the jacobian
    //
    //     The intersection occurs at:
    //       I_y = y_0 - delta * m_y / m_x
    //       I_z = z_0 - delta * 1. / m_x
    //     where: 
    //       delta = xPlane - x_0
    //
    //     The jacobian terms are defined as:
    //        J[i][k] = d I_i / d param_k     where param = { x, y, m_x, m_y }
    //  
    HepMatrix jacobian(2,4);
    //jacobian[0][0] = 0.;                           // dI_y / dx_0  = 0
    //jacobian[0][1] = -delta * inv_m_x;             // dI_y / dm_x  = -d / m_x
    //jacobian[0][2] = 1.;                           // dI_y / dy_0  = 1
    //jacobian[0][3] = delta * m_y * inv_m_x;        // dI_y / dm_y  = -d * m_y / m_x
    
    //jacobian[1][0] = - inv_m_x;                    // dI_z / dx_0  = -1 / m_x
    //jacobian[1][1] = delta * inv_m_x * inv_m_x;    // dI_z / dm_x  = -d / m_x * m_x
    //jacobian[1][2] = 0.;                           // dI_z / dy_0  = 0
    //jacobian[1][3] = 0.;                           // dI_z / dm_y  = 0
    
    jacobian[0][0] = 0.;                           // dI_y / dx_0  = 0
    jacobian[0][1] = 0.;
    jacobian[0][2] = 1.;                           // dI_y / dy_0  = 1
    jacobian[0][3] = 0.;
    
    jacobian[1][0] = - inv_m_x;                    // dI_z / dx_0  = -1 / m_x
    jacobian[1][1] = delta * inv_m_x * inv_m_x;    // dI_z / dm_x  = -d / m_x * m_x
    jacobian[1][2] = 0.;                           // dI_z / dy_0  = 0
    jacobian[1][3] = 0.;                           // dI_z / dm_y  = 0  

    for ( unsigned int i(0);  i < 2; i++ ) {
      for ( unsigned int j(0);  j < 2; j++ ) {
	double currentSum = 0.;
	for ( unsigned int k(1);  k < 5; k++ ) {
	  for ( unsigned int l(1);  l < 5; l++ ) {
	    currentSum += jacobian[i][k-1] * jacobian[j][l-1] * pars(k,l);
	  }
	}
	covAtPlane[i][j] = currentSum;
      }
    }    
  }
  
  void errorAtYPlane(const double delta, const Event::TkrTrackParams& pars, HepMatrix& covAtPlane) {
    
    // get the tk params.  
    // want the x and y and the normal to the plane
    /* not used */ //const double m_x = pars.getxSlope();
    const double m_y = pars.getySlope();
    const double inv_m_y = fabs(m_y) > 1e-9 ? 1./ m_y : 1e9;
    
    // ok input the jacobian
    //
    //     The intersection occurs at:
    //       I_x = x_0 - delta * m_x / m_y
    //       I_z = z_0 - delta * 1. / m_y
    //     where:
    //       delta = yPlane - y_0
    //
    //     The jacobian terms are defined as:
    //        J[i][k] = d I_i / d param_k     where param = { x, y, m_x, m_y }
    //  
    HepMatrix jacobian(2,4);
    //jacobian[0][0] = 1.;                           // dI_x / dx_0  = 1                        
    //jacobian[0][1] = delta * m_x * inv_m_y;        // dI_x / dm_x  = d m_x / m_y
    //jacobian[0][2] = 0.;                           // dI_x / dy_0  = 0
    //jacobian[0][3] = -delta * inv_m_y;             // dI_x / dm_y  = -d / m_y
    
    //jacobian[1][0] = 0.;                           // dI_z / dx_0  = 0                    
    //jacobian[1][1] = 0.;                           // dI_z / dm_x  = 0                    
    //jacobian[1][2] = - inv_m_y;                    // dI_z / dy_0  = -1 / m_y
    //jacobian[1][3] = delta * inv_m_y * inv_m_y;    // dI_z / dm_y  = d / m_y * m_y
    
    jacobian[0][0] = 1.;                           // dI_x / dx_0  = 1                        
    jacobian[0][1] = 0.;
    jacobian[0][2] = 0.;                           // dI_x / dy_0  = 0
    jacobian[0][3] = 0.;
    
    jacobian[1][0] = 0.;                           // dI_z / dx_0  = 0                    
    jacobian[1][1] = 0.;                           // dI_z / dm_x  = 0                    
    jacobian[1][2] = - inv_m_y;                    // dI_z / dy_0  = -1 / m_y
    jacobian[1][3] = delta * inv_m_y * inv_m_y;    // dI_z / dm_y  = d / m_y * m_y  
    
    for ( unsigned int i(0);  i < 2; i++ ) {
      for ( unsigned int j(0);  j < 2; j++ ) {
	double currentSum = 0.;
	for ( unsigned int k(1);  k < 5; k++ ) {
	  for ( unsigned int l(1);  l < 5; l++ ) {
	    currentSum += jacobian[i][k-1] * jacobian[j][l-1] * pars(k,l);
	  }
	}
	covAtPlane[i][j] = currentSum;
      }
    }
  }
  
  void errorAtZPlane(const double /* delta */, const Event::TkrTrackParams& pars, HepMatrix& covAtPlane) {
    
    // ok input the jacobian
    //
    //     The intersection occurs at:
    //       I_x = x_0 - delta * m_x
    //       I_y = y_0 - delta * m_y
    //     where:
    //       delta = zPlane - z_0
    //
    //     The jacobian terms are defined as:
    //        J[i][k] = d I_i / d param_k     where param = { x, y, m_x, m_y }
    //  
    HepMatrix jacobian(2,4);
    
    //jacobian[0][0] = 1.;            // dI_x / dx_0 = 1
    //jacobian[0][1] = -delta;        // dI_x / dm_x = -d
    //jacobian[0][2] = 0.;            // dI_x / dy_0 = 0
    //jacobian[0][3] = 0.;            // dI_x / dm_y = 0
    
    //jacobian[1][0] = 0.;            // dI_y / dx_0 = 0
    //jacobian[1][1] = 0.;            // dI_y / dm_x = 0
    //jacobian[1][2] = 1.;            // dI_y / dy_0 = 1
    //jacobian[1][3] = -delta;        // dI_y / dm_y = -d
    
    jacobian[0][0] = 1.;            // dI_x / dx_0 = 1
    jacobian[0][1] = 0.;        // dI_x / dm_x = -d
    jacobian[0][2] = 0.;            // dI_x / dy_0 = 0
    jacobian[0][3] = 0.;            // dI_x / dm_y = 0
    
    jacobian[1][0] = 0.;            // dI_y / dx_0 = 0
    jacobian[1][1] = 0.;            // dI_y / dm_x = 0
    jacobian[1][2] = 1.;            // dI_y / dy_0 = 1
    jacobian[1][3] = 0.;        // dI_y / dm_y = -d  

    for ( unsigned int i(0);  i < 2; i++ ) {
      for ( unsigned int j(0);  j < 2; j++ ) {
	double currentSum = 0.;
	for ( unsigned int k(1);  k < 5; k++ ) {
	  for ( unsigned int l(1);  l < 5; l++ ) {
	    currentSum += jacobian[i][k-1] * jacobian[j][l-1] * pars(k,l);
	  }
	}
	covAtPlane[i][j] = currentSum;
      }
    }
  }  

  void projectErrorAtPoca(const AcdRecon::TrackData& /* trackData */, const Event::TkrTrackParams& /* trackParams */,
			  const Point& /* poca */, const Vector& /* pocaVector */, double& pocaError) {
    pocaError = 1000.;
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

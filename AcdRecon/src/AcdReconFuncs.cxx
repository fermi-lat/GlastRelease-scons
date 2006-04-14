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

  void tilePlane(const AcdRecon::TrackData& track, const AcdTileDim& tile,
		 double& arcLength, double& localX, double& localY, 
		 double& activeX, double& activeY, double& active2D, HepPoint3D& hitPoint) {

    std::vector<double> dim = tile.dim();
    // Beware: these dimensions are in some sort of local system and for
    // iFace = 1 || 3  x<->y 		
    double dX = dim[0];
    double dY = dim[1];
    double dZ = dim[2];

    const HepPoint3D& center = tile.tileCenter();
    int face = tile.acdId().face();

    crossesPlane(track,center,face,arcLength,localX,localY,hitPoint);

    // If arcLength is negative... had to go backwards to hit plane... 
    if ( arcLength < 0 ) return;
    if(face == 0) {// Top Tile
      activeX = dX/2. - fabs(localX);
      activeY = dY/2. - fabs(localY);
      // Choose which is furthest away from edge (edge @ 0.)
    } else if(face == 1 || face == 3) {// X Side Tile
      activeY = dZ/2  - fabs(localY);
      activeX = dX/2. - fabs(localX);
    } else if(face == 2 || face == 4) {// Y Side Tile
      activeY = dZ/2. - fabs(localY);
      activeX = dX/2. - fabs(localX);
    }
    active2D =  activeX < activeY ? activeX : activeY;
  }

  void tileEdgePoca(const AcdRecon::TrackData& aTrack, const AcdTileDim& tile, 
		    double& arcLength, double& dist, Point& x, Vector& v, int& region) {
    // Loop over all edges allowing only the edges which intersect to form
    // the nearest corner participate.
    // For each pair of corners, make a ray and calculate doca from track
    // Note: this could be done at the end - if no edge solution was found however
    //       doing the full monte does not use much more cpu

    const HepPoint3D* corner = tile.corner();
    dist = 2000.;
    region = -1;
    
    for (int iCorner = 0; iCorner<4; iCorner++) {
      
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

      // Check if x,y,z along edge falls within limits of tile edge.
      double length_2_intersect = raydoca.arcLenRay2();
      if (length_2_intersect > 0 && length_2_intersect < edge_length) {
	double test_dist = raydoca.docaRay1Ray2();
	if ( test_dist < dist ) {
	  dist = test_dist;
	  arcLength = raydoca.arcLenRay1();
	  region = iCorner;
	  x = raydoca.docaPointRay1();
	  v = x - raydoca.docaPointRay2();
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
    const HepPoint3D* corner = tile.corner();
    region = -1;
    dist = 2000;
    
    // First find the nearest corner and the distance to it. 
    double arclen(-1.);    
    double test_dist(2000);
    for (int iCorner = 0; iCorner<4; iCorner++) {      
      HepPoint3D x_isec;
      test_dist = 2000.;
      AcdRecon::pointPoca(aTrack,corner[iCorner],arclen,test_dist,x_isec);
      if ( test_dist < dist ) {
	dist = test_dist;
	region = iCorner+4;
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

  void ribbonPlane(const AcdRecon::TrackData& aTrack, const AcdRibbonDim& ribbon,
		  double& arcLength, double& dist, HepPoint3D& x) {
    
    // get the ribbon id
    const idents::AcdId& acdId = ribbon.acdId();  
    
    // to distiguish X and y ribbons
    static const int ribbonX = 5;

    // this is the value we want to beat
    double best_dist = -2000.;

    // loop over segments
    for ( int segment = 0; segment < 3; segment++ ) {

      // check orientation to determine which segment corresponds to which face
      int face(-1);
      switch ( segment ) {
      case 1: face = 0; break;
      case 0: face = acdId.ribbonOrientation() == ribbonX ? 1 : 2; break;
      case 2: face = acdId.ribbonOrientation() == ribbonX ? 3 : 4; break;
      }

      // Get the beginning and ending points for a line segment
      // that defines a ribbon
      const HepPoint3D& ribbonStartPos = ribbon.ribbonStart()[segment];
      const HepPoint3D& ribbonEndPos = ribbon.ribbonEnd()[segment];
      
      // where does this track cross the plane of the tile 
      HepPoint3D x_isec;
      double test_arc(-1.);
      double localX(0.); double localY(0.);
      AcdRecon::crossesPlane(aTrack,ribbonStartPos,face,test_arc,localX,localY,x_isec);
      if (test_arc < 0) continue;
    
      // Form vector between the beginning of the ribbon and the point where the
      // track intersects the plane of the ribbon
      HepVector3D delta = x_isec - ribbonStartPos;
      // Form a vector for the ribbon
      HepVector3D ribbonVec = ribbonEndPos - ribbonStartPos;
      double prod = delta * ribbonVec.unit();
      // check that the projection of the point to the ribbon occurs within the
      // length of the ribbon segment
      if ((prod < 0) || (prod > ribbonVec.mag())) continue;
    
      double test_dist = sqrt(delta.mag2() - prod*prod);

      // Make this an Active Distance calculation    
      test_dist = ribbon.halfWidth()[segment] - test_dist;

      // test against the best value
      if ( test_dist > best_dist ) {
	dist = test_dist;
	best_dist = test_dist;
	arcLength = test_arc;
	x = x_isec;
      }
    }
  }

  void ribbonPoca(const AcdRecon::TrackData& aTrack, const AcdRibbonDim& ribbon, 
		  double& arcLength, double& dist, Point& x, Vector& v, int& region) {

    // get the ribbon id
    // const idents::AcdId& acdId = ribbon.acdId();  
    
    // this is the value we want to beat
    double best_dist = -2000.;

    // loop over segments
    for ( int segment = 0; segment < 3; segment++ ) {

      // Get the beginning and ending points for a line segment
      // that defines a ribbon
      const HepPoint3D& c1 = ribbon.ribbonStart()[segment];
      const HepPoint3D& c2 = ribbon.ribbonEnd()[segment];
      
      // Will need this to determine limit of the tile edge 
      double edge_length(0.);
      
      // Compute DOCA and DOCA location between the track and edge
      RayDoca raydoca;
      AcdRecon::rayDoca(aTrack,c1,c2,raydoca,edge_length);
      
      // Check if x,y,z along edge falls within limits of tile edge.
      double length_2_intersect = raydoca.arcLenRay2();
      double test_arc = raydoca.arcLenRay1();
      if (length_2_intersect > 0 && test_arc > 0 && length_2_intersect < edge_length) {
	double test_dist = raydoca.docaRay1Ray2();
	// make this an active dist calc
	test_dist = ribbon.halfWidth()[segment] - test_dist;
	if ( test_dist > best_dist ) {
	  dist = test_dist;
	  best_dist = test_dist;
	  arcLength = test_arc;
	  region = segment;
	  x = raydoca.docaPointRay1();
	  v = x - raydoca.docaPointRay2();
	}
      }    
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

  void projectErrorAtPoca(const AcdRecon::TrackData& /* trackData */, const Event::TkrTrackParams& trackParams,
			  const Point& /* poca */, const Vector& /* pocaVector */, double& pocaError) {
    pocaError = 1000.;
  }

}

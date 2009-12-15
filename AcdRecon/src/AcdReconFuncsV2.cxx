#include "../AcdRecon/AcdReconFuncsV2.h"
#include "../AcdRecon/AcdTkrParams.h"

#include "CLHEP/Matrix/Matrix.h"
#include "AcdUtil/AcdTileDim.h"
#include "AcdUtil/AcdRibbonDim.h"
#include "AcdUtil/RayDoca.h"
#include "GlastSvc/Reco/IPropagator.h" 

namespace AcdRecon {

  /**
   * @brief fill the 5 x 4 derivative matrix to get from Tkr covariance form to Acd covariance form
   *
   * @param trackDir the directional cosines
   * @param cov the new covariance matrix
   **/
  void ReconFunctions::fillTkrToAcdCovTranslation(const HepVector3D& dir, HepMatrix& covTrans) {

    covTrans(1,1) = 0.;
    covTrans(1,2) = dir.z() * ( 1. - ( dir.x() * dir.x()) );
    covTrans(1,3) = 0.;
    covTrans(1,4) = -1. * dir.x() * dir.y() * dir.z();
    
    covTrans(2,1) = 0.;
    covTrans(2,2) = -1. * dir.x() * dir.y() * dir.z();
    covTrans(2,3) = 0.;
    covTrans(2,4) = dir.z() * ( 1. - ( dir.y() * dir.y()) );

    covTrans(3,1) = 0.;
    covTrans(3,2) = -1. * dir.x() * dir.z() * dir.z();
    covTrans(3,3) = 0.;
    covTrans(3,4) = -1. * dir.y() * dir.z() * dir.z();

    covTrans(4,1) = 1.;
    covTrans(4,2) = 0.;
    covTrans(4,3) = 0.;
    covTrans(4,4) = 0.;
  
    covTrans(5,1) = 0.;
    covTrans(5,2) = 0.;
    covTrans(5,3) = 1.;
    covTrans(5,4) = 0.;
  }

  /**
   * @brief fill the 4 x 4 covariance matrix from TrkTrackParams
   *
   * @param tkrParams the directional cosines
   * @param cov the  matrix
   **/
  void ReconFunctions::fillCovMatrixFromTkr(const Event::TkrTrackParams& trackParams, HepSymMatrix& cov) {
    for ( int i(1); i < 5; i++) {
      for ( int j(1); j < 5; j++ ){
	cov(i,j) = trackParams(i,j);
      }
    }
  }

  void ReconFunctions::convertToAcdRep(const Event::TkrTrackParams trackParams,
				       double zRef, 
				       AcdRecon::TrackData& acdParams) {

    // first get the parameters
    // Note that Event::TkrTrackParams does indexing from 1.  I think that this is b/c CLHEP hates us.
    
    // Set the reference point
    acdParams.m_point.set( trackParams(Event::TkrTrackParams::xPosIdx), 
			   trackParams(Event::TkrTrackParams::yPosIdx), 
			   zRef );    
    acdParams.m_current = acdParams.m_point;
    
    // get the direction
    // start with the slopes.  Aka tangent vectors
    const double Mx = trackParams(Event::TkrTrackParams::xSlpIdx);
    const double My = trackParams(Event::TkrTrackParams::ySlpIdx);
    
    const double denom2 = 1. + (Mx*Mx) + (My*My);
    const double denom = sqrt(denom2);
    double xDir = Mx/denom;
    double yDir = My/denom;
    double zDir = 1./denom;

    // if this is downgoing flip the direction
    if ( !acdParams.m_upward ) {
      xDir *= -1.;
      yDir *= -1.;
      zDir *= -1.;
    }
    
    acdParams.m_dir.set(xDir,yDir,zDir);

    HepMatrix covTrans(5,4,0);
    fillTkrToAcdCovTranslation( acdParams.m_dir, covTrans );   
    HepSymMatrix tkrCov(4,0);
    fillCovMatrixFromTkr(trackParams, tkrCov);

    // make the covariance matrix in ACD params.
    acdParams.m_cov_orig = tkrCov.similarity(covTrans);
    acdParams.m_cov_prop = acdParams.m_cov_orig;
  }


  bool ReconFunctions::pointPoca(const AcdRecon::TrackData& aTrack, const HepPoint3D& point,
				 double& arcLength, HepPoint3D& poca, HepVector3D& voca ) {
    // Poca to a point.  Done in global frame

    // dX = vector from the Tile Center to the track point
    HepVector3D dX = point - aTrack.m_point;

    // arclength to the poca = dot product between track direction and dX
    arcLength = dX * aTrack.m_dir;

    // don't allow backward solutions 
    if ( arcLength < 0. ) {
      arcLength = 0.01;
    }

    // get the POCA
    poca = aTrack.m_point + arcLength * aTrack.m_dir;   

    // and the VOCA
    voca = poca - point;    
    return true;
  }


  void ReconFunctions::pointPocaError(const AcdRecon::TrackData& track, const HepPoint3D& point,
				      const double& arcLength, const HepVector3D& voca, double& docaError ) {

    // dX = vector from the Tile Center to the current track point
    HepVector3D dX = point - track.m_current;
    
    // build the 3 x 5 derivative matrix to get the cov. on the voca
    HepMatrix B(3,5,0);
    
    B(1,1) = dX.x() * track.m_dir.x() + arcLength;
    B(1,2) = dX.x() * track.m_dir.y();
    B(1,3) = dX.x() * track.m_dir.z();
    B(1,4) = 1. + track.m_dir.x() * track.m_dir.x();
    B(1,5) = track.m_dir.x() * track.m_dir.y();

    B(2,1) = dX.y() * track.m_dir.x();
    B(2,2) = dX.y() * track.m_dir.y() + arcLength;
    B(2,3) = dX.y() * track.m_dir.z();
    B(2,4) = track.m_dir.y() * track.m_dir.x();
    B(2,5) = 1. + track.m_dir.y() * track.m_dir.y();

    B(3,1) = dX.y() * track.m_dir.x();
    B(3,2) = dX.y() * track.m_dir.y();
    B(3,3) = dX.z() * track.m_dir.z() + arcLength;
    B(3,4) = track.m_dir.z() * track.m_dir.x();
    B(3,5) = track.m_dir.z() * track.m_dir.y();

    // The covariance on the VOCA ( use the matrix above )
    HepSymMatrix vocaCov = track.m_cov_prop.similarity( B );

    // The covariance on the DOCA ( project along VOCA )
    HepMatrix V(1,3,0);
    HepVector3D unitVoca = voca.unit();
    V(1,1) = unitVoca.x();
    V(1,2) = unitVoca.y();
    V(1,3) = unitVoca.z();   

    HepSymMatrix docaCov = vocaCov.similarity( V );
    docaError = sqrt( docaCov(1,1) );
  }
  


  bool ReconFunctions::crossesPlane(const AcdRecon::TrackData& track,  
				    const HepPoint3D& planePoint, 
				    const HepVector3D& norm,
				    double& arcLength, HepPoint3D& isec) {
     
    // dX = vector from the current point to the plane center
    HepVector3D dX = planePoint - track.m_point;

    // get the dot products
    double normDotDelta = dX.dot(norm);
    double normDotTrackDir = track.m_dir.dot(norm);

    // Figure out where the on the trajectory this track hits the plane
    arcLength = fabs(normDotTrackDir) > 1e-9 ? (normDotDelta/normDotTrackDir) : -1.;

    // If arcLength is negative... had to go backwards to hit plane... 
    // if ( arcLength < 0 ) return false;

    // Propagate track along by arcLength to get GLOBAL intersection 
    isec = track.m_point + arcLength * track.m_dir;
    return true;
  }

  void ReconFunctions::crossesPlaneError(const AcdRecon::TrackData& track,  
					 const HepPoint3D& planePoint, 
					 const HepTransform3D& toGlobal,
					 const double& arcLength, HepSymMatrix& cov ) {

    // to get the normal vector
    static const HepVector3D zHat(0.,0.,1.);
    HepVector3D norm = toGlobal*zHat;

    // dX = vector from the Tile Center to the current track point
    HepVector3D dX = planePoint - track.m_current;

    // get the dot products
    double normDotDelta = dX.dot(norm);
    double normDotTrackDir = track.m_dir.dot(norm);
    
    double ds_dx1 = -1. * norm.x()/normDotTrackDir;
    double ds_dx2 = -1. * norm.y()/normDotTrackDir;
    double ds_dv1 = -1. * ( normDotDelta * norm.x() ) / ( normDotTrackDir * normDotTrackDir );
    double ds_dv2 = -1. * ( normDotDelta * norm.y() ) / ( normDotTrackDir * normDotTrackDir );
    double ds_dv3 = -1. * ( normDotDelta * norm.z() ) / ( normDotTrackDir * normDotTrackDir );

    // build the 3 x 5 derivative matrix to get the cov. on the voca
    HepMatrix B(3,5,0);
    
    B(1,1) = ds_dv1 * track.m_dir.x() + arcLength;
    B(1,2) = ds_dv2 * track.m_dir.x();
    B(1,3) = ds_dv3 * track.m_dir.x();
    B(1,4) = 1. + ds_dx1 * track.m_dir.x();
    B(1,5) = ds_dx2 * track.m_dir.x();

    B(2,1) = ds_dv1 * track.m_dir.y();
    B(2,2) = ds_dv2 * track.m_dir.y() + arcLength;
    B(2,3) = ds_dv3 * track.m_dir.y();
    B(2,4) = ds_dx1 * track.m_dir.y();
    B(2,5) = 1. + ds_dx2 * track.m_dir.y();

    B(3,1) = ds_dv1 * track.m_dir.z();
    B(3,2) = ds_dv2 * track.m_dir.z();
    B(3,3) = ds_dv3 * track.m_dir.z() + arcLength;
    B(3,4) = ds_dx1 * track.m_dir.z();
    B(3,5) = ds_dx2 * track.m_dir.z();
    
    // The covariance on the intersection point
    HepSymMatrix pointCov = track.m_cov_prop.similarity( B );

    // this is stupid, but hey
    HepMatrix rot(3,3,0);
    rot(1,1) = toGlobal.xx(); rot(1,2) = toGlobal.yx(); rot(1,3) = toGlobal.zx();
    rot(2,1) = toGlobal.xy(); rot(2,2) = toGlobal.yy(); rot(2,3) = toGlobal.zy();
    rot(3,1) = toGlobal.xz(); rot(3,2) = toGlobal.yz(); rot(3,3) = toGlobal.zz();
    
    // The same, in plane ref. frame
    HepSymMatrix planeCov = pointCov.similarity( rot );
    cov(1,1) = planeCov(1,1);    
    cov(1,2) = planeCov(1,2);
    cov(2,2) = planeCov(2,2);

  }


  bool ReconFunctions::rayDoca_withCorner( const AcdRecon::TrackData& track, 
					   const HepPoint3D& c1, const HepPoint3D& c2,
					   double& arcLength, double& rayLength, HepPoint3D& poca, HepVector3D& voca, 
					   int& edgeOrCorner, bool allowCorner ) {

    // Doca between rays (defined by track and two points)  done in global frame
    // Account for finite length of ray defined by two points.
    // If poca occurs beyond end of that ray, returns poca to end point instead
    
    // 
    HepVector3D edge = c2 - c1;
    HepVector3D edgeDir = edge.unit();
    
    // dX = vector from the track ref point to ray ref poit
    HepVector3D dX = track.m_point - c1;

    // Projections of rays along vector between start points
    double d = track.m_dir.dot(dX);
    double e = edgeDir.dot(dX);
    
    // Dot product between two rays to check if parallel
    double b = track.m_dir.dot(edgeDir);
    if ( fabs(b) < 1e-9 ) {
      arcLength = -1;
      return false;
    }
    double f = 1./ ( 1. - b*b);
    arcLength = f* (b*e - d);
    rayLength = f* (e - b*d);
    
    if ( arcLength < 0. ) {
      // Try the value at 0.;
      arcLength = 0.01;
    }

    if ( rayLength < 0 ) {
      // POCA occurs at corner c1
      if ( ! allowCorner ) return false;
      edgeOrCorner = -1;
      pointPoca(track,c1,arcLength,poca,voca);
      return true;
    } else if ( rayLength > edge.mag() ) {      
       // POCA occurs at corner c2
      if ( ! allowCorner ) return false;
      edgeOrCorner = 1;
      pointPoca(track,c2,arcLength,poca,voca);
      return true;
    } else {
      // POCA occurs along ray
      edgeOrCorner = 0;
    }
    
    // Get the Poca and the Voca
    poca = track.m_point + arcLength * track.m_dir;
    HepPoint3D onRay = c1 + rayLength * edgeDir;    
    voca = poca - onRay;
    return true;
  }


  void ReconFunctions::rayDocaError( const AcdRecon::TrackData& track, 
				     const HepPoint3D& c1, const HepPoint3D& c2,
				     const double& arcLength, const double& rayLength, 
				     const HepVector3D& voca, 
				     double& docaError ) {

    // Doca between rays (defined by track and two points)  done in global frame
    // Account for finite length of ray defined by two points.
    // If poca occurs beyond end of that ray, returns poca to end point instead
    
    // 
    HepVector3D edge = c2 - c1;
    HepVector3D edgeDir = edge.unit();
    
    // dX = vector from the current track ref point to ray ref point
    HepVector3D dX = track.m_current - c1;

    // Projections of rays along vector between start points
    double d = track.m_dir.dot(dX);
    double e = edgeDir.dot(dX);
    
    // Dot product between two rays to check if parallel
    double b = track.m_dir.dot(edgeDir);
    double f = 1./ ( 1. - b*b);

    if ( rayLength < 0 ) {
      // POCA occurs at corner c1
      pointPocaError(track,c1,arcLength,voca,docaError);
      return;
    } else if ( rayLength > edge.mag() ) {      
       // POCA occurs at corner c2
      pointPocaError(track,c2,arcLength,voca,docaError);
      return;
    } 

    // get some useful projection vectors
    HepVector3D alpha0 = (e + 2*b*arcLength) * edgeDir - dX;
    HepVector3D alpha1 = (d + 2*b*rayLength) * edgeDir - b*dX;

    HepVector3D beta0 = b*edgeDir - track.m_dir;
    HepVector3D beta1 = edgeDir - b*track.m_dir;

    // build the 3 x 5 derivative matrix to get the cov. on the voca
    HepMatrix B(3,5,0);
        
    B(1,1) = f * ( alpha0.x() * track.m_dir.x() + alpha1.x() * edgeDir.x() ) + arcLength; 
    B(1,2) = f * ( alpha0.x() * track.m_dir.y() + alpha1.x() * edgeDir.y() );
    B(1,3) = f * ( alpha0.x() * track.m_dir.z() + alpha1.x() * edgeDir.z() );
    B(1,4) = 1. + f * ( beta0.x() * track.m_dir.x() + beta0.x() * edgeDir.x() );
    B(1,5) = f * ( beta0.x() * track.m_dir.y() + beta1.x() * edgeDir.y() );

    B(2,1) = f * ( alpha0.y() * track.m_dir.x() + alpha1.y() * edgeDir.x() ); 
    B(2,2) = f * ( alpha0.y() * track.m_dir.y() + alpha1.y() * edgeDir.y() ) + arcLength;
    B(2,3) = f * ( alpha0.y() * track.m_dir.z() + alpha1.y() * edgeDir.z() );
    B(2,4) = f * ( beta0.y() * track.m_dir.x() + beta0.y() * edgeDir.x() );
    B(2,5) = 1. + f * ( beta0.y() * track.m_dir.y() + beta1.y() * edgeDir.y() );

    B(3,1) = f * ( alpha0.z() * track.m_dir.x() + alpha1.z() * edgeDir.x() ); 
    B(3,2) = f * ( alpha0.z() * track.m_dir.y() + alpha1.z() * edgeDir.y() );
    B(3,3) = f * ( alpha0.z() * track.m_dir.z() + alpha1.z() * edgeDir.z() ) + arcLength;
    B(3,4) = f * ( beta0.z() * track.m_dir.x() + beta0.z() * edgeDir.x() );
    B(3,5) = f * ( beta0.z() * track.m_dir.y() + beta1.z() * edgeDir.y() );
   
    // The covariance on the VOCA ( use the matrix above )
    HepSymMatrix vocaCov = track.m_cov_prop.similarity( B );

    // The covariance on the DOCA ( project along voca )
    HepMatrix V(1,3,0);
    HepVector3D unitVoca = voca.unit();
    V(1,1) = unitVoca.x();
    V(1,2) = unitVoca.y();
    V(1,3) = unitVoca.z();   

    HepSymMatrix docaCov = vocaCov.similarity( V );
    docaError = sqrt( docaCov(1,1) );
  }

  void ReconFunctions::tilePlane(const AcdRecon::TrackData& track, const AcdTileDim& tile,
				 AcdRecon::PocaData& data) {
    
    // to get the normal vector
    static const HepVector3D zHat(0.,0.,1.);
    
    // Get all the relevent values for intersections with a tile
    HepPoint3D testHitPoint, testLocal;
    double testActiveX(0.),  testActiveY(0.);

    // loop over volumes of tile
    for (int iVol = 0; iVol < tile.nVol(); iVol++ ) {      
      
      double testArcLength(-1.);
      const AcdTileSection* sect = tile.getSection(iVol);
      HepVector3D norm = sect->m_invTrans * zHat;

      // returns false if intersection is backwards
      if ( ! AcdRecon::ReconFunctions::crossesPlane(track,sect->m_center,norm,
						    testArcLength, testHitPoint) ) continue;

      AcdRecon::ReconFunctions::tilePlaneActiveDistance(tile,iVol,testHitPoint,testLocal,testActiveX,testActiveY);

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
	data.m_volumePlane = iVol;
	data.m_cosTheta = norm.dot(track.m_dir);
      }
    }
    
    int whichVol = data.m_volumePlane;
    if ( whichVol >= 0 ) {
      const AcdTileSection* sect = tile.getSection(whichVol);
      AcdRecon::ReconFunctions::crossesPlaneError(track,sect->m_center,sect->m_invTrans,
						  data.m_arcLengthPlane,data.m_planeError_proj);
      
    }
  }

  
  void ReconFunctions::tilePlaneActiveDistance(const AcdTileDim& tile, int iVol, const HepPoint3D& globalPoint,
					       HepPoint3D& localPoint, double& activeX, double& activeY) {
    const AcdTileSection* sect = tile.getSection(iVol);
    localPoint = sect->m_trans * globalPoint;
    sect->activeDistance(localPoint,iVol,activeX,activeY);
  }


  void ReconFunctions::tileEdgePoca(const AcdRecon::TrackData& track, const AcdTileDim& tile, 
				    AcdRecon::PocaData& data, bool isInside ) {
    // Loop over all edges allowing only the edges which intersect to form
    // the nearest corner participate.
    // For each pair of corners, make a ray and calculate doca from track
    // Note: this could be done at the end - if no edge solution was found however
    //       doing the full monte does not use much more cpu    

    // This is all done in the global frame
    double dist = maxDocaValue;

    // test values
    double testArc(-1.);
    double testRay(-1.);
    HepPoint3D testPoca;
    HepVector3D testVoca;
    int testRegion(-1);

    // loop over volumes of tile
    for (int iVol = 0; iVol < tile.nVol(); iVol++ ) {

      const AcdTileSection* sect = tile.getSection(iVol);
      const HepPoint3D* corner = sect->m_corners;

      for (int iCorner = 0; iCorner<4; iCorner++) {

	// ignored shared edges
	if ( iCorner == sect->m_shared ) continue;

	int ic2 = iCorner == 3 ? 0 : iCorner + 1;
	const HepPoint3D& c1 = corner[iCorner];
	const HepPoint3D& c2 = corner[ic2];  

	// no legal poca found
	if ( ! AcdRecon::ReconFunctions::rayDoca_withCorner(track,c1,c2,
							    testArc,testRay,testPoca,testVoca,
							    testRegion,!isInside) ) continue;
	double testDist = testVoca.mag();
	if ( testDist < dist ) {
	  dist = testDist;
	  // latch values
	  data.m_arcLength = testArc;
	  data.m_rayLength = testRay;
	  data.m_poca = testPoca;
	  data.m_voca = testVoca;
	  data.m_active3D = isInside ? testDist : -1.*testDist;
	  data.m_region = iCorner + (4*(testRegion+1));
	  data.m_volume = iVol;
	}
      }          
    }

    if ( data.m_volume >= 0 ) {
      int iC1 = data.m_region % 4;      
      int iC2 = (data.m_region+1) % 4;

      const AcdTileSection* sect = tile.getSection(data.m_volume);
      const HepPoint3D* corner = sect->m_corners;

      const HepPoint3D& c1 = corner[iC1];
      const HepPoint3D& c2 = corner[iC2]; 
      AcdRecon::ReconFunctions::rayDocaError(track,c1,c2,data.m_arcLength,data.m_rayLength,data.m_voca,
					     data.m_active3DErr_proj);
    }
  }
  

  void ReconFunctions::ribbonPoca(const AcdRecon::TrackData& track, const AcdRibbonDim& ribbon, 
				  AcdRecon::PocaData& data) {

    // 
    double dist = -2000.;
    double testArc(0.);
    double testRay(0.);
    double testDist(0.);
    HepPoint3D testPoca;
    HepVector3D testVoca;
    int testRegion(0);

    int dir(0);
    int startSeg(0);
    int endSeg(0);
    ribbon.getSegmentsIndices(track.m_dir,track.m_upward,startSeg,endSeg,dir);

    for ( int iSeg = startSeg; iSeg != endSeg; iSeg += dir ) {

      const AcdRibbonSegment* seg = ribbon.getSegment(iSeg);
      if ( seg == 0 ) { 
	std::cerr << "No Such Segment " << ribbon.acdId().id()<< ' ' << iSeg << std::endl;
	continue;
      }

      if ( ! AcdRecon::ReconFunctions::rayDoca_withCorner(track,seg->m_start,seg->m_end,
							  testArc,testRay,testPoca,testVoca,testRegion,true) ) continue;

      // Make this an Active Distance calculation  
      testDist = seg->m_halfWidth - testVoca.mag();

      if ( testDist >  dist ) {
	dist = testDist;
	// latch values
	data.m_arcLength = testArc;
	data.m_rayLength = testRay;
	data.m_poca = testPoca;
	data.m_voca = testVoca;
	data.m_active3D = testDist;
	data.m_region = testRegion;
	data.m_volume = iSeg;
      } 
    }

    int whichSeg = data.m_volume;
    if ( whichSeg >= 0 ) {
      data.m_ribbonLength = ribbon.getRibbonLength(whichSeg,data.m_rayLength);
      const AcdRibbonSegment* seg = ribbon.getSegment(whichSeg);
      AcdRecon::ReconFunctions::rayDocaError(track,seg->m_start,seg->m_end,data.m_arcLength,data.m_rayLength,data.m_voca,
					     data.m_active3DErr_proj);
    } else {
      if ( track.m_upward ) {
	std:: cout << "No poca foud " << ribbon.acdId().id() << ' ' 
		   << track.m_point << ' ' << track.m_dir << ' ' << (track.m_upward ? "Up" : "Down" ) << std::endl;
      }
    }
  }

  // initiliaze the kalman propagator
  void ReconFunctions::startPropagator(IPropagator& prop, const Event::TkrTrack& track, const AcdRecon::TrackData& trackData,
				       const double& maxArcLength ) {
    
    // Get the start point & direction of the track & the params & energy also
    const unsigned int hitIndex = trackData.m_upward ? 0 : track.getNumHits() - 1;
    const Event::TkrTrackHit* theHit = track[hitIndex];
    const Point initialPosition = theHit->getPoint(Event::TkrTrackHit::SMOOTHED);
    const Event::TkrTrackParams& trackPars = theHit->getTrackParams(Event::TkrTrackHit::SMOOTHED); 
    
    // setup the propagator
    prop.setStepStart(trackPars,initialPosition.z(),trackData.m_upward); 
    prop.step(maxArcLength);  
  }


  // run the propagtor out to a specified arclength
  void ReconFunctions::propagateToArcLength(IPropagator& prop,
					    const double& arcLength,
					    const AcdRecon::TrackData& trackData,
					    Event::TkrTrackParams& next_params ) {

    Point x_step = prop.getPosition( arcLength );
    next_params = prop.getTrackParams(arcLength,trackData.m_energy,true);
    trackData.m_current.set( x_step.x(),  x_step.y(),  x_step.z() );

    HepMatrix covTrans(5,4,0);
    AcdRecon::ReconFunctions::fillTkrToAcdCovTranslation( trackData.m_dir, covTrans );
    HepSymMatrix tkrCov(4,0);
    AcdRecon::ReconFunctions::fillCovMatrixFromTkr(next_params,tkrCov);
    trackData.m_cov_prop = tkrCov.similarity( covTrans );
  }


  bool ReconFunctions::exitsLat(const AcdRecon::TrackData& trackData,
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

  bool ReconFunctions::entersLat(const AcdRecon::TrackData& trackData, 
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

}

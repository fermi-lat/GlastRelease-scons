#include "AcdPocaTool.h"

#include "AcdUtil/AcdTileDim.h"
#include "AcdUtil/AcdRibbonDim.h"

#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

#include "Event/Recon/AcdRecon/AcdPocaMap.h"
#include "Event/Recon/AcdRecon/AcdHit.h"
#include "Event/Digi/AcdDigi.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GlastSvc/Reco/IPropagator.h" 

#include "CLHEP/Geometry/Transform3D.h"
#include "geometry/Ray.h"
#include "./RayDoca.h"
#include "geometry/Point.h"
#include "geometry/Vector.h"

#include "idents/AcdId.h"

DECLARE_TOOL_FACTORY(AcdPocaTool) ;

const double AcdPocaTool::MaxDoca(2000.);

AcdPocaTool::AcdPocaTool
 ( const std::string & type, 
   const std::string & name,
   const IInterface * parent )
 : AlgTool( type, name, parent )
 { 
   declareInterface<AcdIPocaTool>(this) ; 
   declareProperty("distanceCut",m_distanceCut=200.);
   declareProperty("sigmaCut",m_sigmaCut=5.);
 }

AcdPocaTool::~AcdPocaTool()
{} 

// This function extracts geometry constants from xml
// file using GlastDetSvc
StatusCode AcdPocaTool::initialize()
 {
  MsgStream log(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;
  log<<MSG::INFO<<"BEGIN initialize()"<<endreq ;
 
  return sc;
 }



StatusCode AcdPocaTool::doca (const AcdTileDim& tile,
			      const Event::TkrTrack& aTrack, 
			      PocaData& data) {
  
  // initialize and sanity check
  data.reset(MaxDoca);
  StatusCode sc = tile.statusCode();
  if (sc.isFailure()) return sc;
  const HepPoint3D x0 = aTrack.getInitialPosition();
  const HepVector3D t0 = -(aTrack.getInitialDirection());

  // Get the Tile Center
  const HepPoint3D& xT = tile.tileCenter();
  // dX = vector from the Tile Center to the track start point
  HepVector3D dX = xT - x0;
  // prod = Dot product between dX and the track direction
  double prod = dX * t0;
  // now calculate the doca
  double dist = sqrt(dX.mag2() - prod*prod);
  // compare to the cut value
  if (dist < data.m_dist) {
    // accepted, set the return values 
    data.m_dist = dist;
    data.m_region = 1;
  }
  return sc;
}




StatusCode AcdPocaTool::hitTileDist(const AcdTileDim& tile, 
				    const Event::TkrTrack& aTrack, 
				    PocaData& data) {
  
  // initialize and sanity check
  data.reset(-MaxDoca);  // this is an active dist caculation, want to pick larger numbers
  StatusCode sc = tile.statusCode();
  if (sc.isFailure()) return sc;
  const HepPoint3D x0 = aTrack.getInitialPosition();
  const HepVector3D t0 = -(aTrack.getInitialDirection());
  
  // get the tile geometry
  const idents::AcdId& acdId = tile.acdId();  
  std::vector<double> dim = tile.dim();
  const HepPoint3D& xT = tile.tileCenter();  	
  int iFace = acdId.face();

  // Beware: these dimensions are in some sort of local system and for
  // iFace = 1 || 3  x<->y 		
  double dX = dim[0];
  double dY = dim[1];
  double dZ = dim[2];
  
  // Figure out where in the plane of this face the trajectory hits
  data.m_arcLength = -1.; 
  if(iFace == 0) {// Top Tile. 
    data.m_arcLength = (xT.z()-x0.z())/t0.z();	                
  }
  else if(iFace == 1 || iFace == 3) {// X Side Tile 
    data.m_arcLength = (xT.x()-x0.x())/t0.x();
  }
  else if(iFace == 2 || iFace == 4) {// Y Side Tile
    data.m_arcLength = (xT.y()-x0.y())/t0.y();
  }
  // If arcLength is negative... had to go backwards to hit plane... 
  if( data.m_arcLength < 0.) return sc;
        
  HepPoint3D x_isec = x0 + data.m_arcLength*t0;
  
  HepVector3D local_x0 = x_isec - xT;
  double test_dist(0.);
  if(iFace == 0) {// Top Tile
    double dist_x = dX/2. - fabs(local_x0.x());
    double dist_y = dY/2. - fabs(local_x0.y());	 
    // Choose which is furthest away from edge (edge @ 0.)
    test_dist = (dist_x < dist_y) ? dist_x : dist_y;
  } else if(iFace == 1 || iFace == 3) {// X Side Tile
    double dist_z = dZ/2. - fabs(local_x0.z());
    double dist_y = dX/2. - fabs(local_x0.y());	
    test_dist = (dist_z < dist_y) ? dist_z : dist_y;
  } else if(iFace == 2 || iFace == 4) {// Y Side Tile
    double dist_z = dZ/2. - fabs(local_x0.z());
    double dist_x = dX/2. - fabs(local_x0.x());
    test_dist = (dist_z < dist_x) ? dist_z : dist_x;
  }
  if(test_dist > data.m_dist ) {
    data.m_dist = test_dist;
    data.m_region = 1;
    data.m_poca.set(x_isec.x(),x_isec.y(),x_isec.z());
  }
  return sc;
}

StatusCode AcdPocaTool::tileActiveDist(const AcdTileDim& tile, 
				       const Event::TkrTrack& aTrack, 
				       PocaData& data) {
  
  // initialize and sanity check
  data.reset(-MaxDoca); // this is an active dist caculation, want to pick larger numbers
  StatusCode sc = tile.statusCode();
  if (sc.isFailure()) return sc;
  const HepPoint3D x0 = aTrack.getInitialPosition();
  const HepVector3D t0 = -(aTrack.getInitialDirection());

    // get the tile geometry
  const idents::AcdId& acdId = tile.acdId();  
  std::vector<double> dim = tile.dim();
  const HepPoint3D& xT = tile.tileCenter();  	
  int iFace = acdId.face();

  // Beware: these dimensions are in some sort of local system and for
  // iFace = 1 || 3  x<->y 
  double dX = dim[0];
  double dY = dim[1];
  double dZ = dim[2];
  // The following is required to get the corners to come out correctly
  if(iFace == 1 || iFace == 3) {
    dim[0] = dY;
    dim[1] = dX;
  }
    
  // Figure out where in the plane of this face the trajectory hits
  if(iFace == 0) {// Top Tile. 
    data.m_arcLength = (xT.z()-x0.z())/t0.z();	                
  }
  else if(iFace == 1 || iFace == 3) {// X Side Tile 
    data.m_arcLength = (xT.x()-x0.x())/t0.x();
  }
  else if(iFace == 2 || iFace == 4) {// Y Side Tile
    data.m_arcLength = (xT.y()-x0.y())/t0.y();
  }
  // If arcLength is negative... had to go backwards to hit plane... 
  if(data.m_arcLength < 0.) {    
    data.m_dist = -MaxDoca;
    data.m_region = 0;
    return sc;
  }
  
  HepPoint3D x_isec = x0 + data.m_arcLength*t0;
  
  HepVector3D local_x0 = x_isec - xT;
  if(iFace == 0) {// Top Tile
    double dist_x = dX/2. - fabs(local_x0.x());
    double dist_y = dY/2. - fabs(local_x0.y());	 
    // Choose which is furthest away from edge (edge @ 0.)
    data.m_dist = (dist_x < dist_y) ? dist_x : dist_y;
    data.m_region = (dist_x < dist_y) ? 1 : 2;
  } else if(iFace == 1 || iFace == 3) {// X Side Tile
    double dist_z = dZ/2. - fabs(local_x0.z());
    double dist_y = dX/2. - fabs(local_x0.y());	
    data.m_dist = (dist_z < dist_y) ? dist_z : dist_y;
    data.m_region = (dist_z < dist_y) ? 1 : 2;
  } else if(iFace == 2 || iFace == 4) {// Y Side Tile
    double dist_z = dZ/2. - fabs(local_x0.z());
    double dist_x = dX/2. - fabs(local_x0.x());
    data.m_dist = (dist_z < dist_x) ? dist_z : dist_x;
    data.m_region = (dist_z < dist_x) ? 1 : 2;
  }
  // fall back to the case that this missed the tile
  if (data.m_dist < 0) { 
    sc = docaActiveDist(tile, aTrack, data);
  } else {
    data.m_poca.set(x_isec.x(),x_isec.y(),x_isec.z());
  }
  return sc;
}

StatusCode AcdPocaTool::docaActiveDist(const AcdTileDim& tile,
				       const Event::TkrTrack& aTrack, 
				       PocaData& data) {

  // Initially, we set this to maxDoca, since we are in search of a min value
  // At the end of the method, we negate the final return value
  // initialize and sanity check
  data.reset(MaxDoca);
  StatusCode sc = tile.statusCode();
  if (sc.isFailure()) return sc;
  const HepPoint3D x0 = aTrack.getInitialPosition();
  const HepVector3D t0 = -(aTrack.getInitialDirection());

  // get the tile geometry
  //const idents::AcdId& acdId = tile.acdId();  
  std::vector<double> dim = tile.dim();
  //const HepPoint3D& xT = tile.tileCenter();  	
  //int iFace = acdId.face();

  // Get four corners associated with the tile.
  // Assuming we can avoid, calculation with all 8 corners, 4 should be enough
  // The corners are returned in order (-,-), (-,+), (+,+), (+,-)
  // where third dimension is the one we ignore, since it is associated with
  // tile thickness.
  const HepPoint3D* corner = tile.corner();

  // First find the nearest corner and the distance to it. 
  unsigned int iCorner, i_near_corner = 4; 
  for (iCorner = 0; iCorner<4; iCorner++) {
    HepVector3D dX = corner[iCorner] - x0;
    double prod = dX * t0;
    double dist = sqrt(dX.mag2() - prod*prod);
    if(dist < data.m_dist) {
      data.m_dist = dist;
      i_near_corner = iCorner;
      data.m_region = 3;
      data.m_arcLength = prod;
      HepPoint3D x_isec = x0 + prod*t0;
      data.m_poca.set(x_isec.x(),x_isec.y(),x_isec.z());
      HepVector3D pocaVect = x_isec - corner[iCorner];
      data.m_pocaVector.set(pocaVect.x(),pocaVect.y(),pocaVect.z());
    }
  }

  // Loop over all edges allowing only the edges which intersect to form
  // the nearest corner participate.
  // For each pair of corners, make a ray and calculate doca from track
  // Note: this could be done at the end - if no edge solution was found however
  //       doing the full monte does not use much more cpu
  Point trackPos(x0.x(), x0.y(), x0.z());
  Vector trackDir(t0.x(), t0.y(), t0.z());
  Ray track(trackPos, trackDir);
  for (iCorner = 0; iCorner<4; iCorner++) {
    
    //  WBA: Naively I thought one could limit the edges investigated - not so!
    //	if(iCorner != i_near_corner && (iCorner+1)%4 != i_near_corner) continue;
    
    Point pos0(corner[iCorner].x(), corner[iCorner].y(), 
	       corner[iCorner].z());
    Point pos1;
    if(iCorner==3)
      pos1.set(corner[0].x(), corner[0].y(), corner[0].z());
    else
      pos1.set(corner[iCorner+1].x(), corner[iCorner+1].y(), 
	       corner[iCorner+1].z());
    Vector dir = pos1 - pos0;
    Ray edge(pos0, dir);
    
    // Will need this to determine limit of the tile edge 
    double edge_length = dir.magnitude();
    
    // Compute DOCA and DOCA location between the track and edge
    RayDoca raydoca = RayDoca(track, edge);

    // Check if x,y,z along edge falls within limits of tile edge.
    double length_2_intersect = raydoca.arcLenRay2();
    if (length_2_intersect > 0 && length_2_intersect < edge_length) {
      double test_dist = raydoca.docaRay1Ray2();
      if ( test_dist < data.m_dist ) {
	data.m_dist = test_dist;
	data.m_arcLength = raydoca.arcLenRay1();
	data.m_region = 1;
	data.m_poca = raydoca.docaPointRay1();
	data.m_pocaVector = data.m_poca - raydoca.docaPointRay2();
      }
    }    
  }

  // Negate the return distance, because we call this function in the case
  // where the track fails to pierce a tile
  data.m_dist *= -1;

  return sc;

}

StatusCode AcdPocaTool::hitRibbonDist(const AcdRibbonDim& ribbon,
				      const Event::TkrTrack& aTrack, 
				      PocaData& data) {
  //Purpose and Method:  Calculate ActiveDistance for Ribbons
  // To simplify the situation - we assume the ribbons are merely lines 
  // Since the ribbons are segmented in the geometry - we treat the side segments and 
  // top segments separately, taking the minimum distance over all hit segmen

  // initialize and sanity check
  data.reset(-MaxDoca);  // this is an active dist caculation, want to pick larger numbers
  StatusCode sc = ribbon.statusCode();
  if (sc.isFailure()) return sc;
  const HepPoint3D x0 = aTrack.getInitialPosition();
  const HepVector3D t0 = -(aTrack.getInitialDirection());

  // get the ribbon id
  const idents::AcdId& acdId = ribbon.acdId();  

  // Need to reconstruct the required volume ids to retrieve the geometry information
  // For now we brute force it.. and construct what we need
  int topSegment;
  int sideFace[2];

  static const int ribbonX = 5;

  // check orientation to determine which segment number to retrieve for top
  if ( acdId.ribbonOrientation() == ribbonX) {
    topSegment = 1;
    // ribbons that are along x-axis on the top go down faces 1,3
    sideFace[0] = 1; sideFace[1] = 3;
  } else {
    topSegment = 2;
    // ribbons that are along the y-axis on the top go down faces 2,4
    sideFace[0] = 2; sideFace[1] = 4;
  }
  
  // Now loop over 3 segments for this ribbon
  // We want to grab the 2 side segments and retrieve their dimensions
  // Also want to grap one of the top segments and for the orientation that is cut
  // into 5 segments - we want to take the middle segment and use that to construct
  // a line covering the whole top - we'll use that line to calculate active distance
  int isegment;
  for (isegment = 0; isegment < 3; isegment++) {
 
    // After all that, we have the beginning and ending points for a line segment
    // that defines a ribbon
    const HepPoint3D& ribbonStartPos = ribbon.ribbonStart()[isegment];
    const HepPoint3D& ribbonEndPos = ribbon.ribbonEnd()[isegment];
    
    int iFace = isegment == 1 ? 0 : sideFace[isegment/2]; 
    
    // Figure out where in the plane of this face the trajectory hits
    double arc_dist = -1.; 
    if(iFace == 0) {// Top ribbon
      arc_dist = (ribbonStartPos.z()-x0.z())/t0.z();	                
    } else if(iFace == 1 || iFace == 3) {// X Side ribbon
      arc_dist = (ribbonStartPos.x()-x0.x())/t0.x();
    } else if(iFace == 2 || iFace == 4) {// Y Side ribbon
      arc_dist = (ribbonStartPos.y()-x0.y())/t0.y();
    }
    
    if (arc_dist < 0) continue;
    
    // position of hit in the plane of the ribbon
    HepPoint3D x_isec = x0 + arc_dist*t0;
    
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
    test_dist = ribbon.halfWidth()[isegment] - test_dist;
    if( test_dist > data.m_dist ) {
      data.m_dist = test_dist;
      data.m_arcLength = arc_dist;
      data.m_region = 1;
      data.m_poca.set(x_isec.x(),x_isec.y(),x_isec.z());
    } // endif
  }
  return sc;
}


StatusCode AcdPocaTool::makePoca(const Event::TkrTrack& aTrack, int iTrack,
				 const PocaData& poca, const idents::AcdId& acdId,
				 IPropagator& g4PropTool, Event::AcdTkrPoca*& thePoca) {


  StatusCode sc = StatusCode::SUCCESS;
  thePoca = 0;
  if ( poca.m_dist < -1.*m_distanceCut ) return sc;

  Event::TkrTrackParams paramsAtPoca;
  sc = getParamsAtPoca(aTrack,true,poca.m_arcLength,g4PropTool,paramsAtPoca);
  if ( sc.isFailure() ) {
    return sc;
  }

  double pocaError(1000.);
  sc = projectError(poca,paramsAtPoca,pocaError);
  if ( sc.isFailure() ) {
    return sc;
  }
  
  double sigma = poca.m_dist / pocaError;
  if ( sigma < -1.*m_sigmaCut ) return sc;
  
  thePoca = new Event::AcdTkrPoca(acdId,iTrack,poca.m_dist,pocaError,poca.m_region,poca.m_poca,paramsAtPoca);
  return sc;
					      
}

StatusCode AcdPocaTool::getParamsAtPoca(const Event::TkrTrack& aTrack, bool forward, double arcLength,
					IPropagator& g4PropTool,
					Event::TkrTrackParams& paramsAtPoca) {

  StatusCode sc = StatusCode::SUCCESS;
  
  // get the first (or last) hit
  const unsigned int hitIndex = forward ? 0 : aTrack.getNumHits() - 1;  
  const Event::TkrTrackHit* theHit = aTrack[hitIndex];
  if ( theHit == 0 ) return StatusCode::FAILURE;

  // initialize the propagator
  const Point initialPosition = theHit->getPoint(Event::TkrTrackHit::SMOOTHED);
  const Event::TkrTrackParams& trackPars = theHit->getTrackParams(Event::TkrTrackHit::SMOOTHED);
  g4PropTool.setStepStart(trackPars,initialPosition.z(),forward);   
  double startEnergy = theHit->getEnergy();

  // get the params at the poca
  paramsAtPoca = g4PropTool.getTrackParams(arcLength,startEnergy,true);
  return sc;
 
}

StatusCode AcdPocaTool::projectError(const PocaData& /* poca */, const Event::TkrTrackParams& /* paramsAtPoca */,
				     double& pocaError) {

  StatusCode sc = StatusCode::SUCCESS;
  pocaError = 10000.;

  return sc;
  
}

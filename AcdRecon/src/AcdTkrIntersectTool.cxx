#include "AcdTkrIntersectTool.h"

#include "AcdIPocaTool.h"
#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

#include "Event/Recon/TkrRecon/TkrTrackHit.h"

#include "CLHEP/Matrix/Matrix.h"

#include "geometry/Point.h"
#include "geometry/Vector.h"

#include "../AcdRecon/AcdReconFuncs.h"

#include "idents/AcdId.h"

#include "CLHEP/Geometry/Transform3D.h"

// Define the fiducial volume of the LAT
// FIXME -- this should come for some xml reading service
//
// top is defined by planes at + 754.6 -> up to stacking of tiles
// sides are defined by planes at +-840.14
// the bottom of the ACD is at the z=-50 plane

// Later we add 10 cm to make sure that we catch everything
const double AcdTkrIntersectTool::s_top_distance  = 754.6;     // center of tiles in cols 1 and 3
const double AcdTkrIntersectTool::s_side_distance  = 840.14;   // center of tiles in sides
const double AcdTkrIntersectTool::s_bottom_distance  = -50.;   // bottom of ACD

DECLARE_TOOL_FACTORY(AcdTkrIntersectTool)

AcdTkrIntersectTool::AcdTkrIntersectTool
 ( const std::string & type, 
   const std::string & name,
   const IInterface * parent )
 : AlgTool( type, name, parent )
 { declareInterface<AcdITkrIntersectTool>(this) ; }

AcdTkrIntersectTool::~AcdTkrIntersectTool()
{} 

// This function extracts geometry constants from xml
// file using GlastDetSvc
StatusCode AcdTkrIntersectTool::initialize()
 {
  MsgStream log(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;
  log<<MSG::INFO<<"BEGIN initialize()"<<endreq ;

  // GlastDetSvc
  sc = service("GlastDetSvc", m_detSvc);
  if (sc.isFailure()) {
    log<< MSG::ERROR<<"GlastDetSvc not found"<<endreq ;
    return StatusCode::FAILURE ;
  } 
  
  // G4Propagator
  if( ! toolSvc()->retrieveTool("G4PropagationTool", m_G4PropTool)) {
    log << MSG::ERROR << "Couldn't find the G4Propagator!" << endreq;
    return StatusCode::FAILURE;
  }

  // AcdPocaTool  
  if( ! toolSvc()->retrieveTool("AcdPocaTool", m_pocaTool)) {
    log << MSG::ERROR << "Couldn't find the AcdPocaTool!" << endreq;
    return StatusCode::FAILURE;
  }

  log<<MSG::INFO<<"END initialize()"<<endreq ;
  return StatusCode::SUCCESS ;
 }

StatusCode  AcdTkrIntersectTool::makeIntersections(IPropagator& prop, 
						   const AcdRecon::TrackData& track,
						   const AcdRecon::ExitData& data,						   
						   const AcdRecon::PocaDataPtrMap& pocaMap,
						   const AcdRecon::AcdHitMap& hitMap,
						   AcdGeomMap& geomMap,
						   Event::AcdTkrIntersectionCol& intersections,
						   Event::AcdTkrGapPocaCol& gapPocas) {
  
  MsgStream log(msgSvc(),name()) ;
  
  // figure out the direction
  bool forward = data.m_arcLength > 0.;
  
  // setup some of the variables for the step loop
  idents::VolumeIdentifier volId;
  idents::VolumeIdentifier prefix = m_detSvc->getIDPrefix();

  // storage inside the loop
  double arcLength(0.);

  // latch which intersection we use for the gap
  AcdRecon::PocaData* ribbonDataForGap(0);
  AcdRecon::PocaData* tileDataForGap(0);
  double bestTileActDist(-2000.);

  // this is in case we have to make any new pocas
  std::list<AcdRecon::PocaData> ownedPocaData;

  // loop over the GEANT steps
  int numSteps = m_G4PropTool->getNumberSteps();  
  int istep(0);
  for(istep = 0; istep < numSteps; ++istep) {

    // running sums of arcLength, this has to be before the check on volume id
    double pathLength = m_G4PropTool->getStepArcLen(istep); 
    arcLength  += pathLength;

    // check to see if this step is in the Acd
    volId = prop.getStepVolumeId(istep);
    volId.prepend(prefix);
    if ( ! checkVolId(volId) ) continue;
    
    // which face of detector (0=top,1=-x,2=-y,3=+x,4=+y)
    int face = volId[1];
    idents::AcdId acdId = idents::AcdId(volId);
 
    // propgate the postion and error matrix to the intersection
    Point x_step = prop.getPosition( arcLength );
    Event::TkrTrackParams next_params = prop.getTrackParams(arcLength,track.m_energy,true);
    
    // query the detector service for the local frame    
    HepTransform3D transform; 
    HepPoint3D center(0., 0., 0.);
    StatusCode sc = m_detSvc->getTransform3DByID(volId, &transform);
      
    if ( sc.isFailure() ) {
      MsgStream log(msgSvc(),name()) ;
      log << MSG::ERROR << "Could not get transform for id " << volId.name() << endreq;
      return sc;
    }        
    HepPoint3D xT = transform * center;

    // was this tile hit or not
    AcdRecon::AcdHitMap::const_iterator findHit = hitMap.find(acdId);
    unsigned char hitMask = findHit != hitMap.end() ? (unsigned char)(findHit->second) : 0; 

    AcdRecon::PocaData* pocaData(0);
    AcdRecon::PocaDataPtrMap::const_iterator itrFind = pocaMap.find(acdId);
    if ( itrFind != pocaMap.end() ) {
      // already did this tile, just get the values we already computed
      pocaData = itrFind->second;
    } else {            
      // haven't done this tile yet.  
      ownedPocaData.push_back( AcdRecon::PocaData() );
      pocaData = &( ownedPocaData.back() );
      pocaData->m_id = acdId;
      if ( acdId.ribbon() ) { 	
	const AcdRibbonDim* ribbon = geomMap.getRibbon(acdId,*m_detSvc);
	if ( ribbon->statusCode().isFailure() ) {
	  log << MSG::ERROR << "Failed to get geom for a ribbon " << acdId.id() << ' ' << volId.name() << endreq;
	  return StatusCode::FAILURE;
	}
	HepPoint3D p;
	AcdRecon::ribbonPlane(track,*ribbon,
			      pocaData->m_arcLengthPlane,pocaData->m_active2D,p);
	pocaData->m_inPlane.set(p.x(),p.y(),p.z());
	AcdRecon::ribbonPoca(track,*ribbon,
			     pocaData->m_arcLength,pocaData->m_active3D,pocaData->m_poca,pocaData->m_pocaVector,pocaData->m_region);
      } else if ( acdId.tile() ) {
	const AcdTileDim* tile = geomMap.getTile(acdId,*m_detSvc);
	if ( tile->statusCode().isFailure() ) {
	  log << MSG::ERROR << "Failed to get geom for a tile " << acdId.id() << ' ' << volId.name() << endreq;
	  return StatusCode::FAILURE;
	}
	HepPoint3D p;
	AcdRecon::tilePlane(track,*tile,
			    pocaData->m_arcLengthPlane,pocaData->m_localX,pocaData->m_localY,
			    pocaData->m_activeX,pocaData->m_activeY,pocaData->m_active2D,p);
	pocaData->m_inPlane.set(p.x(),p.y(),p.z());
	AcdRecon::tileEdgePoca(track,*tile,
			       pocaData->m_arcLength,pocaData->m_active3D,pocaData->m_poca,pocaData->m_pocaVector,pocaData->m_region);
	sc = holePoca(track,*pocaData,*tile,gapPocas);
	if ( sc.isFailure() ) return sc;
      }
    }
    AcdRecon::projectToPlane(track,next_params,face,xT,*pocaData);

    if ( acdId.ribbon() ) {
      // cull out the top Y - ribbons
      if ( face != 0 || acdId.ribbonOrientation() != 5 ) {
	ribbonDataForGap = pocaData;
      }
    } else if ( acdId.tile() ) {
      float testDist = face == 0 ? pocaData->m_activeY : pocaData->m_activeX;
      if ( testDist > bestTileActDist ) {
	tileDataForGap = pocaData;
	bestTileActDist = testDist;
      }    
    }
        
    HepMatrix localCovMatrix(2,2);
    localCovMatrix[0][0] = pocaData->m_localCovXX;
    localCovMatrix[1][1] = pocaData->m_localCovYY;
    localCovMatrix[0][1] = localCovMatrix[1][0] = pocaData->m_localCovXY;

    double localPosition[2]; 
    localPosition[0] = pocaData->m_localX;
    localPosition[1] = pocaData->m_localY;
    
    // ok, we have everything we need, build the AcdTkrIntersection 
    Event::AcdTkrIntersection* iSect = 
      new Event::AcdTkrIntersection(acdId,track.m_index,
				    x_step,
				    localPosition,localCovMatrix,
				    (forward ? arcLength : -1.*arcLength), pathLength,
				    hitMask, pocaData->m_cosTheta);
    
    // print to the log if debug level is set
    iSect->writeOut(log);
    
    // append to collection and increment counter
    intersections.push_back(iSect);
  }

  if ( ribbonDataForGap != 0 ) {
    gapPocaRibbon(track,data,*ribbonDataForGap,gapPocas);
  } else if ( tileDataForGap != 0 ) {
    gapPocaTile(track,data,*tileDataForGap,gapPocas);
  } else {
    fallbackToNominal(track,data,gapPocas);
  } 

  return StatusCode::SUCCESS ;
}


StatusCode AcdTkrIntersectTool::exitsLAT(const Event::TkrTrack& aTrack, bool upward,
					 AcdRecon::ExitData& data) {
 
  // which end of track to get hit from
  const unsigned int hitIndex = upward ? 0 : aTrack.getNumHits() - 1;
  const Event::TkrTrackHit* theHit = aTrack[hitIndex];

  // get position, direction
  const Point initialPosition = theHit->getPoint(Event::TkrTrackHit::SMOOTHED);
  const Vector initialDirection = upward ? 
    -1.* theHit->getDirection(Event::TkrTrackHit::SMOOTHED) : 
    theHit->getDirection(Event::TkrTrackHit::SMOOTHED);
  
  return exitsLAT(initialPosition,initialDirection,upward,data);
}
 

StatusCode AcdTkrIntersectTool::exitsLAT(const Point& initialPosition, const Vector& initialDirection, bool upward,
					 AcdRecon::ExitData& data) {

  MsgStream log(msgSvc(),name()) ;   

  // hits -x or +x side ?
  const double normToXIntersection =  initialDirection.x() > 0 ?  
    -1.*s_side_distance - initialPosition.x() :    // hits -x side
    1.*s_side_distance - initialPosition.x();      // hits +x side  
  const double slopeToXIntersection = fabs(initialDirection.x()) > 1e-9 ? 
    -1. / initialDirection.x() : (normToXIntersection > 0. ? 1e9 : -1e9);
  const double sToXIntersection = normToXIntersection * slopeToXIntersection;
  
  // hits -y or +y side ?
  const double normToYIntersection = initialDirection.y() > 0 ?  
    -1.*s_side_distance - initialPosition.y() :    // hits -y side
    1.*s_side_distance - initialPosition.y();      // hits +y side
  const double slopeToYIntersection = fabs(initialDirection.y()) > 1e-9 ? 
    -1. / initialDirection.y() : (normToYIntersection > 0. ? 1e9 : -1e9); 
  const double sToYIntersection = normToYIntersection * slopeToYIntersection;
  
  // hits top or bottom
  const double normToZIntersection = upward ? 
    s_top_distance - initialPosition.z() :
    s_bottom_distance - initialPosition.z();
  const double slopeToZIntersection = 1. / initialDirection.z();
  const double sToZIntersection = normToZIntersection * slopeToZIntersection;
  
  if ( (upward && initialDirection.z() < 0) || (!upward && initialDirection.z() > 0) ) {
    log << MSG::ERROR << "Downgoing track " << upward << ' ' << initialDirection.z()  << endreq;
  }

  // pick the closest place
  if ( sToXIntersection < sToYIntersection ) {
    if ( sToXIntersection < sToZIntersection ) {
      // hits X side
      data.m_arcLength = sToXIntersection;
      data.m_face = initialDirection.x() > 0 ? 1 : 3;
    } else {
      // hits Z side
      data.m_arcLength = sToZIntersection;
      data.m_face = upward ? 0 : 5;
    }     
  } else {
    if ( sToYIntersection < sToZIntersection ) {
      // hits X side
      data.m_arcLength = sToYIntersection;
      data.m_face = initialDirection.y() > 0 ? 2 : 4;
    } else {
      // hits Z side
      data.m_arcLength = sToZIntersection;
      data.m_face = upward ? 0 : 5;
    }     
  }

  // protect against negative arcLengths
  if ( data.m_arcLength < 0. ) {
    log << MSG::ERROR << "Negative Arclength to intersection " << data.m_arcLength << endreq;
    return StatusCode::FAILURE;
  }

  // extrapolate to the i-sect
  data.m_x = initialPosition - data.m_arcLength*initialDirection;

  // flip the sign of the arclength for downgoing side
  data.m_arcLength *= upward ? 1. : -1.;

  return StatusCode::SUCCESS;
}
 
bool AcdTkrIntersectTool::checkVolId(idents::VolumeIdentifier& volId) {
  if ( ! volId.isAcd() ) return false;
  switch ( volId[2] ) {
  case 40: // Acd Tile
  case 41: // Acd ribbon
    if ( volId.size() < 5 ) {
      // FIXME -- need to verify that this is really dead material
      return false;
    }      
    break;
  default:
    // Dead Acd material
    return false;
  }    
  return true;
}

StatusCode AcdTkrIntersectTool::fallbackToNominal(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
						  Event::AcdTkrGapPocaCol& gapPocas) {
  // FIXME -- do this right  
  static std::vector<float> gapListX;
  static std::vector<float> gapListY;
  static bool ready(false);

  if ( !ready ) {
    gapListX.push_back(-835.);    gapListY.push_back(-835.);
    gapListX.push_back(-502.);    gapListY.push_back(-502.);
    gapListX.push_back(-166.);    gapListY.push_back(-168.);
    gapListX.push_back(166.);    gapListY.push_back(168.);
    gapListX.push_back(502.);    gapListY.push_back(502.);
    gapListX.push_back(835.);    gapListY.push_back(835.);
    ready = true;
  }

  double hit(0.);
  std::vector<float>* gapList(0);

  // Make a pocaData struct to pass the info around
  AcdRecon::PocaData pocaData;

  AcdRecon::AcdGapType gapType(AcdRecon::None);

  switch ( data.m_face ) {
  case 0:
    // top
    hit = data.m_x.x();
    gapList = &gapListX;
    break;
  case 1:
  case 3:
    hit = data.m_x.y();
    gapList = &gapListY;
    break;
  case 2:
  case 4:
    hit = data.m_x.x();
    gapList = &gapListX;
    break;
  case 5:
    return StatusCode::SUCCESS ;
    break;
  }

  // compare to the known gaps
  float maxDoca(-2000.);
  int whichGap(-1);
  
  for ( unsigned int i(0); i < gapList->size(); i++ ) {
    float doca = 5. - fabs( hit - (*gapList)[i] );
    if ( doca > maxDoca ) {
      maxDoca = doca;
      whichGap = i;
    }
  }

  switch ( data.m_face ) {
  case 0:
    // top
    gapType = ( whichGap == 0 || whichGap == 5 ) ? AcdRecon::TopCornerEdge : AcdRecon::Y_RibbonTop;
    break;
  case 1:
  case 3:
    gapType = ( whichGap == 0 || whichGap == 5 ) ? AcdRecon::SideCornerEdge : AcdRecon::Y_RibbonSide;
    break;
  case 2:
  case 4:
    gapType = ( whichGap == 0 || whichGap == 5 ) ? AcdRecon::SideCornerEdge : AcdRecon::X_RibbonSide;
    break;
  default:
    break;
  }
  int row(0); int col(0);
  int gapIndex = whichGap -1;
  if ( gapIndex == -1 ) gapIndex = 0;
  if ( gapIndex == 4 ) gapIndex = 1;


  idents::AcdGapId gapId(gapType,gapIndex,data.m_face,row,col);

  if ( false ) {
    std::cout << "ToGap: [" << data.m_x.x() << ' ' <<  data.m_x.y() << ' ' <<  data.m_x.z() << "] " 
	      << gapType << ' ' << data.m_face << ' ' << gapIndex << ' ' << maxDoca << std::endl;
  }

  Event::AcdTkrGapPoca* poca(0);
  StatusCode sc = makeGapPoca(gapId,track,pocaData,maxDoca,poca);
  if ( sc.isFailure() ) return sc;

  gapPocas.push_back(poca);
  return sc;
}
  
StatusCode AcdTkrIntersectTool::holePoca(const AcdRecon::TrackData& /* track */, 
					 const AcdRecon::PocaData& /* pocaData */, const AcdTileDim& /* tile */,
					 Event::AcdTkrGapPocaCol& /* gapPocas */) {

  return StatusCode::SUCCESS ;
}

StatusCode AcdTkrIntersectTool::gapPocaRibbon(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
					      const AcdRecon::PocaData& pocaData, Event::AcdTkrGapPocaCol& gapPocas) {

  AcdRecon::AcdGapType gapType(AcdRecon::None);  
  unsigned char face = (unsigned char)(data.m_face);
  unsigned char gap = 0;
  switch ( face ) {
  case 0:
    gapType = AcdRecon::Y_RibbonTop;
    gap = 2;
    break;
  case 1:
  case 3:
    gapType = AcdRecon::Y_RibbonSide;
    gap = 2;
    break;
  case 2:
  case 4:
    gapType = AcdRecon::X_RibbonSide;
    gap = 2;
    break;
  case 5:
    return StatusCode::SUCCESS;
  }

  unsigned char col = (unsigned char) ( pocaData.m_id.ribbonNum() );
  unsigned char row = 0;
  
  idents::AcdGapId gapId(gapType,gap,face,row,col);

  double distance = pocaData.m_active2D;

  Event::AcdTkrGapPoca* poca(0);
  StatusCode sc = makeGapPoca(gapId,track,pocaData,distance,poca);
  if ( sc.isFailure() ) return sc;

  gapPocas.push_back(poca);
  return StatusCode::SUCCESS;
  
}
  
StatusCode AcdTkrIntersectTool::gapPocaTile(const AcdRecon::TrackData& track, const AcdRecon::ExitData& /* data */,
					    const AcdRecon::PocaData& pocaData, Event::AcdTkrGapPocaCol& gapPocas) {
  
  signed char gapType = (unsigned char)(AcdRecon::None);
  unsigned char face = (unsigned char)(pocaData.m_id.face());
  unsigned char row = (unsigned char)(pocaData.m_id.row());
  unsigned char col = (unsigned char)(pocaData.m_id.column());
  unsigned char gap = 0;

  double distance(0.);

  // bottow side tile have gaps only at the ends
  if ( face != 0 && row == 3 ) {
    gapType = AcdRecon::SideCornerEdge;
    if ( pocaData.m_localX > 0 ) { gap = 1; }    
    distance = -1.*pocaData.m_activeX;
  } else {
    if ( pocaData.m_localX > 0 ) { col++; }
    if ( pocaData.m_localY > 0 ) { row++; }    
    switch ( face ) {
    case 0:
      gapType = ( col == 0 || col == 5 ) ?  AcdRecon::TopCornerEdge : AcdRecon::Y_RibbonTop;    
      if ( pocaData.m_localY > 0 ) gap = 1;
      distance = -1.*pocaData.m_activeX;
      break;
    case 1:
    case 3:    
      gapType = ( col == 0 || col == 5 ) ? AcdRecon::SideCornerEdge : AcdRecon::Y_RibbonSide;
      if ( pocaData.m_localX > 0 ) gap = 1;
      distance = -1.*pocaData.m_activeX;
      break;
    case 2:
    case 4:
      gapType = ( col == 0 || col == 5 ) ? AcdRecon::SideCornerEdge : AcdRecon::X_RibbonSide;
      if ( pocaData.m_localX > 0 ) gap = 1;
      distance = -1.*pocaData.m_activeX;
      break;
    case 5:
      return StatusCode::SUCCESS;
    }
    if ( pocaData.m_localX > 0 ) { col--; }
    if ( pocaData.m_localY > 0 ) { row--; }
  }

 
  idents::AcdGapId gapId(gapType,gap,face,row,col);
 
  Event::AcdTkrGapPoca* poca(0);
  StatusCode sc = makeGapPoca(gapId,track,pocaData,distance,poca);
  if ( sc.isFailure() ) return sc;

  gapPocas.push_back(poca);
  return StatusCode::SUCCESS;
}

StatusCode AcdTkrIntersectTool::makeGapPoca(idents::AcdGapId& gapId, const AcdRecon::TrackData& track, const AcdRecon::PocaData& pocaData,
					    double distance, Event::AcdTkrGapPoca*& poca) {

  poca = 0;
  double arcLength = track.m_upward ? pocaData.m_arcLength : -1* pocaData.m_arcLength;
    
  float local[2];
  local[0] = pocaData.m_activeX;
  local[1] = pocaData.m_activeY;
  HepMatrix localCov(2,2);
  localCov[0][0] = pocaData.m_localCovXX;
  localCov[1][1] = pocaData.m_localCovYY;
  localCov[0][1] = localCov[1][0] = pocaData.m_localCovXY;

  // temp storage
  static Event::AcdTkrLocalCoords localCoords;
  static Event::AcdPocaData pd;

  localCoords.set(local,pocaData.m_path,pocaData.m_cosTheta,pocaData.m_region,localCov);
  pd.set(arcLength,distance,pocaData.m_active3DErr,pocaData.m_poca,pocaData.m_pocaVector);
  poca = new Event::AcdTkrGapPoca(gapId,track.m_index,localCoords,pd);

  return StatusCode::SUCCESS;
}

StatusCode AcdTkrIntersectTool::makeTkrPoint(const AcdRecon::TrackData& /*track*/, const AcdRecon::ExitData& data,
					     const Event::TkrTrackParams& params, Event::AcdTkrPoint*& tkrPoint ) {
  tkrPoint = 0;
  
  tkrPoint = new Event::AcdTkrPoint(/* track.m_index, */ data.m_arcLength,data.m_face,data.m_x,params);
  return StatusCode::SUCCESS;
}


StatusCode AcdTkrIntersectTool::entersLAT(const Point& initialPosition, const Vector& initialDirection, bool upward,
					  AcdRecon::ExitData& data) {

  MsgStream log(msgSvc(),name()) ;   

  // where does the track start relative to +-X sides
  // and how long before it hits one of the sides
  // this evals to -1 if between sides
  double sToXIntersection(-1.);
  if ( fabs(initialPosition.x()) > s_side_distance ) {
    double normToXIntersection = initialPosition.x() < 0 ? 
      -1.*s_side_distance - initialPosition.x() :    // hits -x side first
      1.*s_side_distance - initialPosition.x();      // hits +x side frist
    const double slopeToXIntersection = fabs(initialDirection.x()) > 1e-9 ? 
      1. / initialDirection.x() : (normToXIntersection > 0. ? 1e9 : -1e9);
    sToXIntersection = normToXIntersection * slopeToXIntersection;
    // propagate to that point, make sure that other two values inside LAT also
    if ( sToXIntersection > 0 ) {
      Point xPlaneInter = initialPosition;  xPlaneInter += sToXIntersection* initialDirection;
      if (  fabs(xPlaneInter.y()) < s_side_distance  &&
	    xPlaneInter.z() > s_bottom_distance &&
	    xPlaneInter.z() < s_top_distance ) {
	data.m_arcLength = sToXIntersection;
	data.m_x = xPlaneInter;
	data.m_face = initialPosition.x() < 0 ? 1 : 3;
      }
    }
  }  

  // where does the track start relative to +-Y sides
  // this evals to -1 if between sides
  double sToYIntersection(-1.);
  if ( fabs(initialPosition.y()) > s_side_distance ) {
    double normToYIntersection = initialPosition.y() < 0 ? 
      -1.*s_side_distance - initialPosition.y() :    // hits -y side first
      1.*s_side_distance - initialPosition.y();      // hits +y side frist
    const double slopeToYIntersection = fabs(initialDirection.y()) > 1e-9 ? 
      1. / initialDirection.y() : (normToYIntersection > 0. ? 1e9 : -1e9);
    sToYIntersection = normToYIntersection * slopeToYIntersection;
    if ( sToYIntersection > 0 && 
	 ( sToYIntersection < data.m_arcLength || data.m_arcLength < 0 ) ) {    
      // propagate to that point, make sure that other two values inside LAT also
      Point yPlaneInter = initialPosition;  yPlaneInter += sToYIntersection* initialDirection;
      if (  fabs(yPlaneInter.x()) < s_side_distance  &&
	    yPlaneInter.z() > s_bottom_distance &&
	    yPlaneInter.z() < s_top_distance ) {
	data.m_arcLength = sToYIntersection;
	data.m_x = yPlaneInter;
	data.m_face = initialPosition.y() < 0 ? 2 : 4;	
      }
    }
  }  

  // where does the track start relative to +-Z sides
  // this evals to -1 if between sides
  double sToZIntersection(-1.);
  if ( initialPosition.z() < s_bottom_distance ||
       initialPosition.z() > s_top_distance ) {
    double normToZIntersection = initialPosition.z() < 0 ? 
      s_bottom_distance - initialPosition.z() :    // hits -z side first
      s_top_distance - initialPosition.z();      // hits +z side frist
    const double slopeToZIntersection = fabs(initialDirection.z()) > 1e-9 ? 
      1. / initialDirection.z() : (normToZIntersection > 0. ? 1e9 : -1e9);
    sToZIntersection = normToZIntersection * slopeToZIntersection;
    if ( sToZIntersection > 0 && 
	 ( sToZIntersection < data.m_arcLength || data.m_arcLength < 0 ) ) {    
      // propagate to that point, make sure that other two values inside LAT also
      Point zPlaneInter = initialPosition;  zPlaneInter += sToZIntersection* initialDirection;
      if (  fabs(zPlaneInter.x()) < s_side_distance  &&
	    fabs(zPlaneInter.y()) < s_side_distance ) {
	data.m_arcLength = sToZIntersection;
	data.m_x = zPlaneInter;
	data.m_face = initialPosition.z() > s_top_distance ? 0 : 5;	
      }
    }
  }  

  return StatusCode::SUCCESS;
}

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
#include "../AcdRecon/AcdTkrParams.h"
#include "./AcdPocaSorter.h"

#include "idents/AcdId.h"

#include "CLHEP/Geometry/Transform3D.h"

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
  
  // GlastDetSvc
  sc = service("AcdGeometrySvc", m_acdGeomSvc);
  if (sc.isFailure()) {
    log<< MSG::ERROR<<"AcdGeometrySvc not found"<<endreq ;
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

  // Sort out the list of existing POCAs
  AcdPocaSorter pocaSorter(track.m_upward? AcdPocaSorter::Upward: AcdPocaSorter::Downward, pocaMap);

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
  Event::TkrTrackParams next_params;
  AcdTkrParams paramsAtISect;
  std::vector<AcdPocaSorter::AcdPocaHolder> pocasCrossed;

  for(istep = 0; istep < numSteps; ++istep) {

    // running sums of arcLength, this has to be before the check on volume id
    double pathLength = m_G4PropTool->getStepArcLen(istep); 
    arcLength  += pathLength;

    // check to see if this step is in the Acd
    volId = prop.getStepVolumeId(istep);

    AcdRecon::PocaData* pocaData(0);
    unsigned nPocaCrossed(0);
    unsigned char hitMask(0);

    idents::AcdId acdId;
    int face(-1);

    if ( checkVolId(volId) ) {
      
      volId.prepend(prefix);
      acdId = idents::AcdId(volId);
      face = volId[1];

      // get the ordered list of POCA up to this arclength
      // This resets the pocasCrossed list before filling
      nPocaCrossed = pocaSorter.getPocasToArclength(arcLength,pocasCrossed);    
      StatusCode sc = fillPocaData(track,volId,pocaMap,geomMap,ownedPocaData,pocaData);
      if ( sc.isFailure() ) {
	log << MSG::ERROR << "AcdTkrIntersectionTool::fillPocaData failed" << endreq;
	return sc;
      }
      pocasCrossed.push_back(AcdPocaSorter::AcdPocaHolder(AcdPocaSorter::PlanePoca,pocaData));
      nPocaCrossed++;
      // query the detector service for the local frame    
      // was this tile hit or not
      AcdRecon::AcdHitMap::const_iterator findHit = hitMap.find(acdId);
      hitMask = findHit != hitMap.end() ? (unsigned char)(findHit->second) : 0; 
    } else {
      if ( istep == numSteps-1 ) {
	nPocaCrossed = pocaSorter.getPocasToArclength(arcLength,pocasCrossed);    
      }
    }

    if ( nPocaCrossed == 0 ) continue;

    for ( std::vector<AcdPocaSorter::AcdPocaHolder>::iterator itrPoca = pocasCrossed.begin(); itrPoca != pocasCrossed.end(); itrPoca++ ) {
      Point x_step = prop.getPosition( itrPoca->arclength() );
      next_params = prop.getTrackParams(arcLength,track.m_energy,true);

      AcdRecon::PocaData* nextPoca = itrPoca->poca();

      // protect against poorly estimated volumes
      if ( nextPoca->m_volume < 0 ) continue;

      if ( nextPoca->m_id.ribbon() ) {
	;
      } else if ( nextPoca->m_id.tile() ) {	

	const AcdTileDim* tile = geomMap.getTile(nextPoca->m_id,*m_acdGeomSvc);
	// propgate the postion and error matrix to the intersection
	AcdRecon::propagateToArcLength(prop,track,arcLength,next_params,paramsAtISect);	
	AcdRecon::projectErrorToPlane(paramsAtISect,tile->localFrameVectors(nextPoca->m_volume),nextPoca->m_planeError);
	AcdRecon::projectErrorToPocaVector(paramsAtISect,nextPoca->m_pocaVector,nextPoca->m_active3DErr);	
      }    
    }
    
    if ( pocaData == 0 ) continue;

    // Back to the intersection at hand
    if ( acdId.ribbon() ) {
      // cull out the top Y - ribbons
      if ( face != 0 || acdId.ribbonOrientation() != 5 ) {
	ribbonDataForGap = pocaData;
      }
    } else if ( acdId.tile() ) {
      const AcdTileDim* tile = geomMap.getTile(acdId,*m_acdGeomSvc);
      if ( holePoca(track,*pocaData,*tile,gapPocas).isFailure() ) {
	return StatusCode::FAILURE;
      }

      // Check to see if this is the Tile to use for getting the distance to the 
      // closest gap
      float testDist = face == 0 ? pocaData->m_activeY : pocaData->m_activeX;
      if ( testDist > bestTileActDist ) {
	tileDataForGap = pocaData;
	bestTileActDist = testDist;
      }    
    }

    double localPosition[2]; 
    localPosition[0] = pocaData->m_inPlane.x();
    localPosition[1] = pocaData->m_inPlane.y();
    
    // ok, we have everything we need, build the AcdTkrIntersection 
    Event::AcdTkrIntersection* iSect = 
      new Event::AcdTkrIntersection(acdId,track.m_index,
				    pocaData->m_hitsPlane,
				    localPosition,pocaData->m_planeError,
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
    gapPocaTile(track,data,*tileDataForGap,geomMap,gapPocas);
  } else {
    // Didn't hit any tiles or ribbons
    // Check to see if it exited bottom (data.m_face==5)
    // Otherwise computed distance to nominal gap
    if ( data.m_face != 5 ) fallbackToNominal(track,data,gapPocas);
  } 
  
  gapPocaCorner(track,data,gapPocas);

  return StatusCode::SUCCESS ;
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
	      << gapType << ' ' << data.m_face << ' ' << gapIndex << ' ' << maxDoca << ' ' << data.m_arcLength << std::endl;
  }

  Event::AcdTkrGapPoca* poca(0);
  StatusCode sc = makeGapPoca(gapId,track,pocaData,maxDoca,poca);
  if ( sc.isFailure() ) return sc;

  gapPocas.push_back(poca);
  return sc;
}
  

StatusCode AcdTkrIntersectTool::fillPocaData(const AcdRecon::TrackData& track, const idents::VolumeIdentifier& volId,
					     const AcdRecon::PocaDataPtrMap& pocaMap, AcdGeomMap& geomMap,
					     std::list<AcdRecon::PocaData>& ownedPocaData, AcdRecon::PocaData*& pocaData) {  

  // Which ACD id is this?
  idents::AcdId acdId = idents::AcdId(volId);
  /* int face = volId[1]; */
  
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
      const AcdRibbonDim* ribbon = geomMap.getRibbon(acdId,*m_acdGeomSvc);
      if ( ribbon->statusCode().isFailure() ) {
	MsgStream log(msgSvc(),name()) ;
	log << MSG::ERROR << "Ribbon statusCode is failure" << endreq;
	return StatusCode::FAILURE;
      }
      HepPoint3D p;
      AcdRecon::ribbonPlane(track,*ribbon,*pocaData);
      AcdRecon::ribbonPoca(track,*ribbon,
			     pocaData->m_arcLength,pocaData->m_activeY,
			     pocaData->m_active3D,pocaData->m_poca,pocaData->m_pocaVector,pocaData->m_region);
      } else if ( acdId.tile() ) {
	const AcdTileDim* tile = geomMap.getTile(acdId,*m_acdGeomSvc);
	if ( tile->statusCode().isFailure() ) {
	  MsgStream log(msgSvc(),name()) ;
	  log << MSG::ERROR << "Tile statusCode is failure" << endreq;
	  return StatusCode::FAILURE;
	}
	HepPoint3D p;
	AcdRecon::tilePlane(track,*tile,*pocaData);
	AcdRecon::tileEdgePoca(track,*tile,
			       pocaData->m_arcLength,pocaData->m_active3D,pocaData->m_poca,pocaData->m_pocaVector,
			       pocaData->m_region);
	
      }
    }    

  return StatusCode::SUCCESS;
}

StatusCode AcdTkrIntersectTool::holePoca(const AcdRecon::TrackData& /* track */, 
					 const AcdRecon::PocaData& /* pocaData */, const AcdTileDim& /* tile */,
					 Event::AcdTkrGapPocaCol& /* gapPocas */) {

  return StatusCode::SUCCESS ;
}

StatusCode AcdTkrIntersectTool::gapPocaRibbon(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
					      const AcdRecon::PocaData& pocaData, Event::AcdTkrGapPocaCol& gapPocas) {

  AcdRecon::AcdGapType gapType(AcdRecon::None);  
  int face = data.m_face;
  unsigned gap = 0;
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

  unsigned col = pocaData.m_id.ribbonNum();
  unsigned row = 0;
  
  idents::AcdGapId gapId(gapType,gap,face,row,col);

  double distance = pocaData.m_active3D;

  Event::AcdTkrGapPoca* poca(0);
  StatusCode sc = makeGapPoca(gapId,track,pocaData,distance,poca);

  if ( sc.isFailure() ) return sc;

  gapPocas.push_back(poca);
  return StatusCode::SUCCESS;
  
}
  
StatusCode AcdTkrIntersectTool::gapPocaTile(const AcdRecon::TrackData& track, const AcdRecon::ExitData& /* data */,
					    const AcdRecon::PocaData& pocaData, AcdGeomMap& geomMap, Event::AcdTkrGapPocaCol& gapPocas) {
  
  signed char gapType = (unsigned char)(AcdRecon::None);
  unsigned char face = (unsigned char)(pocaData.m_id.face());
  unsigned char row = (unsigned char)(pocaData.m_id.row());
  unsigned char col = (unsigned char)(pocaData.m_id.column());
  unsigned char gap = 0;

  double distance(0.);

  AcdRecon::PocaData ribbonPocaData;

  // bottow side tile have gaps only at the ends
  if ( face != 0 && row == 3 ) {
    gapType = AcdRecon::SideCornerEdge;
    if ( pocaData.m_inPlane.x() > 0 ) { gap = 1; }    
    distance = -1.*pocaData.m_activeX;
  } else {
    if ( pocaData.m_inPlane.y() > 0 ) { row++; }    
    switch ( face ) {
    case 0: // TOP, Normal col and row counting 
      if ( pocaData.m_inPlane.x() > 0 ) { col++; }
      gapType = ( col == 0 || col == 5 ) ?  AcdRecon::TopCornerEdge : AcdRecon::Y_RibbonTop;    
      if ( pocaData.m_inPlane.y() > 0 ) gap = 1;
      distance = -1.*pocaData.m_activeX;
      break;
    case 1:  // -Y face, column go opposite of local x
      if ( pocaData.m_inPlane.x() < 0 ) { col++; gap = 1; }
      gapType = ( col == 0 || col == 5 ) ? AcdRecon::SideCornerEdge : AcdRecon::X_RibbonSide;      
      distance = -1.*pocaData.m_activeX;
      break;
    case 3:  // +Y face, column go with local x
      if ( pocaData.m_inPlane.x() > 0 ) { col++; gap = 1; }
      gapType = ( col == 0 || col == 5 ) ? AcdRecon::SideCornerEdge : AcdRecon::X_RibbonSide;
      distance = -1.*pocaData.m_activeX;
      break;
    case 2:  // -X face, column go opposite of local x
      if ( pocaData.m_inPlane.x() < 0 ) { col++; gap = 1; }
      gapType = ( col == 0 || col == 5 ) ? AcdRecon::SideCornerEdge : AcdRecon::Y_RibbonSide;
      distance = -1.*pocaData.m_activeX;
      break;
    case 4:  // +X face, column go with local x
      if ( pocaData.m_inPlane.x() > 0 ) { col++; gap = 1; }
      gapType = ( col == 0 || col == 5 ) ? AcdRecon::SideCornerEdge : AcdRecon::Y_RibbonSide;
      distance = -1.*pocaData.m_activeX;
      break;
    case 5:
      return StatusCode::SUCCESS;
    }

    idents::AcdId whichRibbon;
    if ( gapType == AcdRecon::Y_RibbonTop || gapType == AcdRecon::Y_RibbonSide ) {
      whichRibbon = idents::AcdId(6,col-1);
    } else if ( gapType == AcdRecon::X_RibbonSide ) {
      whichRibbon = idents::AcdId(5,col-1);      
    }
    if ( pocaData.m_inPlane.x() > 0 ) { 
      switch ( face ) {
      case 0:
      case 3:
      case 4:
	col--;
	break;
      case 1:
      case 2:
	col++;
	break;
      }
    }
    if ( pocaData.m_inPlane.y() > 0 ) { row--; }

    if ( whichRibbon.ribbon() ) {
      const AcdRibbonDim* ribbon = geomMap.getRibbon(whichRibbon,*m_acdGeomSvc);      
      AcdRecon::ribbonPoca(track,*ribbon,ribbonPocaData.m_arcLength,ribbonPocaData.m_ribbonLength,
			   ribbonPocaData.m_active3D,ribbonPocaData.m_poca,ribbonPocaData.m_pocaVector,ribbonPocaData.m_region);
    }
  }
   
  idents::AcdGapId gapId(gapType,gap,face,row,col);

  Event::AcdTkrGapPoca* poca(0);
  StatusCode sc = makeGapPoca(gapId,track,ribbonPocaData,distance,poca);
  if ( sc.isFailure() ) return sc;

  gapPocas.push_back(poca);
  return StatusCode::SUCCESS;
}

StatusCode AcdTkrIntersectTool::gapPocaCorner(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
					      Event::AcdTkrGapPocaCol& gapPocas) {

  // iterate over all corner gaps
  double bestDist(2000.);
  int whichCorner(-1);

  AcdRecon::PocaData pocaData;

  unsigned int iCorner(0);
  for (iCorner=0; iCorner<4; iCorner++) {
    const Ray& gapRay = m_acdGeomSvc->getCornerGapRay(iCorner);
    // Compute DOCA between the track and gap ray 

    double testArcLen(-1.);
    double testRayLen(-1.);
    double testDist(-2000.);
    Point testPoint;
    Vector testDir;
    int testRegion(-1);

    AcdRecon::rayDoca_withCorner(track,gapRay,testArcLen,testRayLen,testDist,testPoint,testDir,testRegion);

    // only take forward intersections
    if ( testArcLen < 0. ) continue;
    if ( testDist < bestDist ) {
      bestDist = testDist;
      pocaData.m_arcLength = testArcLen;
      float sign = ( (track.m_point.x() * track.m_dir.y()) - (track.m_point.y() * track.m_dir.x()) ) > 0 ? 1. : -1;
      pocaData.m_active3D = testDist * sign;
      pocaData.m_poca = testPoint;
      pocaData.m_pocaVector = testDir;
      pocaData.m_region = testRegion;
      whichCorner = iCorner;
    }
  }

  idents::AcdGapId gapId(AcdRecon::CornerRay,whichCorner+4,data.m_face,0,0);
  
  Event::AcdTkrGapPoca* poca(0);
  StatusCode sc = makeGapPoca(gapId,track,pocaData,pocaData.m_active3D,poca);
  if ( sc.isFailure() ) return sc;

  gapPocas.push_back(poca);
  return StatusCode::SUCCESS;
}

StatusCode AcdTkrIntersectTool::makeGapPoca(idents::AcdGapId& gapId, const AcdRecon::TrackData& track, const AcdRecon::PocaData& pocaData,
					    double distance, Event::AcdTkrGapPoca*& poca) {

  poca = 0;
  double arcLength = track.m_upward ? pocaData.m_arcLength : -1* pocaData.m_arcLength;
    
  float local[2];
  switch ( gapId.gapType() ) {
  case AcdRecon::X_RibbonSide:
  case AcdRecon::Y_RibbonSide: 
  case AcdRecon::Y_RibbonTop: 
    // this gap has been calculated from a ribbon
    local[0] = pocaData.m_active3D;
    local[1] = pocaData.m_ribbonLength;
    break;
  default:
    // this gap has been calculated from a tile or a corner ray
    local[0] = pocaData.m_activeX;
    local[1] = pocaData.m_activeY;
  }

  // temp storage
  static Event::AcdTkrLocalCoords localCoords;
  static Event::AcdPocaData pd;

  localCoords.set(local,pocaData.m_path,pocaData.m_cosTheta,pocaData.m_region,pocaData.m_planeError);
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


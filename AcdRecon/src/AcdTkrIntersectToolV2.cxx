#include "AcdTkrIntersectToolV2.h"

#include "AcdIPocaTool.h"
#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

#include "Event/Recon/TkrRecon/TkrTrackHit.h"

#include "CLHEP/Matrix/Matrix.h"

#include "geometry/Point.h"
#include "geometry/Vector.h"

#include "../AcdRecon/AcdReconFuncsV2.h"
#include "./AcdPocaSorter.h"

#include "idents/AcdId.h"

#include "CLHEP/Geometry/Transform3D.h"

#include <TMath.h>

DECLARE_TOOL_FACTORY(AcdTkrIntersectToolV2)

AcdTkrIntersectToolV2::AcdTkrIntersectToolV2
 ( const std::string & type, 
   const std::string & name,
   const IInterface * parent )
   : AlgTool( type, name, parent ),
     m_pocaTool(0),
     m_G4PropTool(0),
     m_detSvc(0),
     m_acdGeomSvc(0)
 { declareInterface<AcdITkrIntersectToolV2>(this) ; 
 }

AcdTkrIntersectToolV2::~AcdTkrIntersectToolV2()
{} 

// This function extracts geometry constants from xml
// file using GlastDetSvc
StatusCode AcdTkrIntersectToolV2::initialize()
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
  if( ! toolSvc()->retrieveTool("AcdPocaToolV2", m_pocaTool)) {
    log << MSG::ERROR << "Couldn't find the AcdPocaToolV2!" << endreq;
    return StatusCode::FAILURE;
  }

  log<<MSG::INFO<<"END initialize()"<<endreq ;
  return StatusCode::SUCCESS ;
 }

StatusCode  AcdTkrIntersectToolV2::makeIntersections(IPropagator& prop, 
						     const Event::TkrTrackParams& trackParams,
						     const AcdRecon::TrackData& track,
						     const AcdRecon::ExitData& data,
						     const AcdRecon::AcdHitMap& hitMap,
						     AcdGeomMap& geomMap,
						     AcdRecon::PocaDataPtrMap& pocaMap,
						     std::list<AcdRecon::PocaData> ownedPocaData,
						     std::vector<Event::AcdTkrGapPoca*>& gapPocas) {
  
  MsgStream log(msgSvc(),name()) ;
  
  // figure out the direction
  bool forward = data.m_arcLength > 0.;
  
  // first figure out how far to extrapolate track
  double maxArcLength = forward ? data.m_arcLength : -1. * data.m_arcLength;
  
  // setup some of the variables for the step loop
  idents::VolumeIdentifier volId;
  idents::VolumeIdentifier prefix = m_detSvc->getIDPrefix();

  // Sort out the list of existing POCAs
  AcdPocaSorter pocaSorter(track.m_upward? AcdPocaSorter::Upward: AcdPocaSorter::Downward, pocaMap);

  double largestArc = pocaSorter.finalArc();
  if ( largestArc > maxArcLength ) {
    maxArcLength = largestArc;
  }
  maxArcLength += 20.;  
  
  // protect against negative arcLengths
  if ( forward && maxArcLength < 0 ) {
    log << MSG::ERROR << "Negative Arclength to upper intersection " << maxArcLength << endreq;
    return StatusCode::FAILURE;
  }
  if ( !forward && maxArcLength < 0 ) {
    log << MSG::ERROR << "Negative Arclength to lower intersection " << maxArcLength << endreq;
    return StatusCode::FAILURE;
  }
  
  AcdRecon::ReconFunctions::startPropagator(prop,trackParams,track,maxArcLength);

  // storage inside the loop
  double arcLength(0.);

  // latch which intersection we use for the gap
  AcdRecon::PocaData* tileDataForGap(0);
  double bestTileActDist(-2000.);

  // this is in case we have to make any new pocas

  // loop over the GEANT steps
  int numSteps = m_G4PropTool->getNumberSteps();  
  int istep(0);
  Event::TkrTrackParams next_params;
  std::vector<AcdPocaSorter::AcdPocaHolder> pocasCrossed;

  unsigned totalPoca(0);

  for(istep = 0; istep < numSteps; ++istep) {

    // running sums of arcLength, this has to be before the check on volume id
    double pathLength = m_G4PropTool->getStepArcLen(istep); 
    arcLength  += pathLength;

    // check to see if this step is in the Acd
    volId = prop.getStepVolumeId(istep);

    // This is the data associated with the step across the struck ACD Element
    // It is only filled for steps across elements
    AcdRecon::PocaData* pocaDataCrossed(0);

    unsigned nPocaCrossed(0);
    unsigned char hitMask(0);

    idents::AcdId acdId;

    if ( checkVolId(volId) ) {
      
      volId.prepend(prefix);
      acdId = idents::AcdId(volId);

      // get the ordered list of POCA up to this arclength
      // This resets the pocasCrossed list before filling
      nPocaCrossed = pocaSorter.getPocasToArclength(arcLength,pocasCrossed);   

      // go ahead and fill the poca for this acdID
      StatusCode sc = fillPocaData(track,volId,pocaMap,arcLength,geomMap,ownedPocaData,pocaDataCrossed);
      if ( sc.isFailure() ) {	
	log << MSG::ERROR << "AcdTkrIntersectionTool::fillPocaData failed to get poca data" << std::endl
	    << data.m_x.x() << ' ' << data.m_x.y() << ' ' << data.m_x.z() << ' ' << data.m_arcLength << endreq;
	pocaDataCrossed = 0;
      } else {
	// add a new poca holder object, just of completeness
	pocasCrossed.push_back(AcdPocaSorter::AcdPocaHolder(AcdPocaSorter::PlanePoca,pocaDataCrossed));
	nPocaCrossed++;
      }
      // was this tile hit or not
      AcdRecon::AcdHitMap::const_iterator findHit = hitMap.find(acdId);      
      hitMask = findHit != hitMap.end() ? getHitMask(*(findHit->second)) : 0;
    } else {
      // last step. grab all the pocas we have done yet
      if ( istep >= numSteps-1 ) {
	nPocaCrossed = pocaSorter.getPocasToArclength(maxArcLength,pocasCrossed);    
      }
    }
    
    if ( nPocaCrossed == 0 ) { continue; }
    totalPoca += nPocaCrossed;

    for ( std::vector<AcdPocaSorter::AcdPocaHolder>::iterator itrPoca = pocasCrossed.begin(); 
	  itrPoca != pocasCrossed.end(); itrPoca++ ) {

      AcdRecon::PocaData* nextPoca = itrPoca->poca();
      // propgate the postion and error matrix to the intersection
      AcdRecon::ReconFunctions::propagateToArcLength(prop,itrPoca->arclength(),track,next_params);
      if ( nextPoca->m_id.ribbon() ) {
	int ivol = nextPoca->m_volume;
	if ( ivol < 0 ) { 
	  log << MSG::ERROR << "Bad ribbon volume " << nextPoca->m_id.id() << ' ' << nextPoca->m_volume << endreq;	  
	  continue;
	}
	const AcdRibbonDim* ribbon = geomMap.getRibbon(nextPoca->m_id,*m_acdGeomSvc);
	const AcdRibbonSegment* seg = ribbon->getSegment(ivol);
	// get the relative arc length, 
	double relArcLength = 0.;
	if ( itrPoca->type() == AcdPocaSorter::PlanePoca ) {
	  HepPoint3D cent = seg->m_start + ( 0.5 * seg->m_vect );	  
	  AcdRecon::ReconFunctions::crossesPlaneError(track,cent,seg->m_invTrans,
				      relArcLength,nextPoca->m_planeError_prop);
	} else {
	  const HepPoint3D& c1 = seg->m_start;
	  const HepPoint3D& c2 = seg->m_end;	  
	  AcdRecon::ReconFunctions::rayDocaError(track,c1,c2,
						 relArcLength,nextPoca->m_rayLength,
						 nextPoca->m_voca,nextPoca->m_active3DErr_prop);
	}
      } else if ( nextPoca->m_id.tile() ) {	
	// protect against poorly estimated volumes
	int ivol = nextPoca->m_volume;
	if ( ivol < 0 ) { 
	  log << MSG::ERROR << "Bad tile volume " << nextPoca->m_id.id() << ' ' << nextPoca->m_volume << endreq;	  
	  continue;
	}
	const AcdTileDim* tile = geomMap.getTile(nextPoca->m_id,*m_acdGeomSvc);
	const AcdTileSection* sect = tile->getSection(ivol);
	
	// Check to see if this is the Tile to use for getting the distance to the 
	// closest gap
	float testDist = nextPoca->m_activeX;
	
	if ( testDist > bestTileActDist && nextPoca->m_activeY > 0 ) {
	  tileDataForGap = nextPoca;
	  bestTileActDist = testDist;
	}    

	// get the relative arc length
	double relArcLength = 0.;
	if ( itrPoca->type() == AcdPocaSorter::PlanePoca ) {
	  AcdRecon::ReconFunctions::crossesPlaneError(track,sect->m_center,sect->m_invTrans,
						      relArcLength,nextPoca->m_planeError_prop);
	} else {
	  const HepPoint3D& c1 = sect->m_corners[nextPoca->m_region % 4];
	  const HepPoint3D& c2 = sect->m_corners[(nextPoca->m_region+1) % 4];	  
	  AcdRecon::ReconFunctions::rayDocaError(track,c1,c2,
						 relArcLength,nextPoca->m_rayLength,
						 nextPoca->m_voca,nextPoca->m_active3DErr_prop);
	}
      }    
    }

    // check to see if this was a step across an ACD element, if not, go on to next step
    if ( pocaDataCrossed == 0 ) continue;

    if ( acdId.tile() ) {
      AcdTileDim* tile = const_cast<AcdTileDim*>(geomMap.getTile(acdId,*m_acdGeomSvc));
      if ( holePoca(track,*pocaDataCrossed,*tile,gapPocas).isFailure() ) {
	return StatusCode::FAILURE;
      }
    }
  }

  if ( ! pocaSorter.isDone() ) {
    log <<  MSG::ERROR << "AcdPocaSorter did not finish " << (forward ? 'f' : 'b') << ' ' 
	<< (pocaSorter.dir() == AcdPocaSorter::Upward ? 'u' : 'd') << ' ' 
	<< arcLength << ' ' << pocaSorter.finalArc() << ' ' << maxArcLength << ' ' 
	<< totalPoca << ' ' << pocaSorter.size() << endreq;
    pocaSorter.runOut();
  }

  AcdRecon::ReconFunctions::propagateToArcLength(prop,arcLength,track,next_params);

  if ( tileDataForGap != 0 ) {
    gapPocaTile(track,data,*tileDataForGap,geomMap,gapPocas);
  } else {
    // Didn't hit any tiles or ribbons
    // Check to see if it exited bottom (data.m_face==5)
    // Otherwise computed distance to nominal gap
    if ( data.m_face != 5 && data.m_x.z() > 0 ) {
      log << MSG::ERROR << "Didn't find any tile intersections for track that does not exit bottom of ACD" << std::endl 
	  << "Exit Point: " << data.m_face << ' ' << data.m_x.x() << ' ' << data.m_x.y() << ' ' << data.m_x.z() <<  endreq;
    }
  }   

  return StatusCode::SUCCESS ;
}
 

bool AcdTkrIntersectToolV2::checkVolId(idents::VolumeIdentifier& volId) {
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

unsigned char AcdTkrIntersectToolV2::getHitMask(const Event::AcdHit& aHit) {
  unsigned char hitMask = 0;
  hitMask |= aHit.getAcceptMapBit(Event::AcdHit::A) ? 0x1 : 0;
  hitMask |= aHit.getAcceptMapBit(Event::AcdHit::B) ? 0x2 : 0;
  hitMask |= aHit.getHitMapBit(Event::AcdHit::A) ? 0x4 : 0;
  hitMask |= aHit.getHitMapBit(Event::AcdHit::B) ? 0x8 : 0;
  return hitMask;
}

StatusCode AcdTkrIntersectToolV2::fillPocaData(const AcdRecon::TrackData& track, const idents::VolumeIdentifier& volId,
					       const AcdRecon::PocaDataPtrMap& pocaMap, const double& arc, 
					       AcdGeomMap& geomMap,
					       std::list<AcdRecon::PocaData>& ownedPocaData, 
					       AcdRecon::PocaData*& pocaData) {  
  
  // Which ACD id is this?
  idents::AcdId acdId = idents::AcdId(volId);
  AcdRecon::PocaDataPtrMap::const_iterator itrFind = pocaMap.find(acdId);
  if ( itrFind != pocaMap.end() ) {
    // already did this tile, just get the values we already computed
    pocaData = itrFind->second;    
    return StatusCode::SUCCESS;
  }
   
  // haven't done this tile yet.  
  HepPoint3D thePoint = track.m_point + (arc *  track.m_dir);

  MsgStream log(msgSvc(),name()) ;
  log << MSG::WARNING << "Missed Poca for ID " << acdId.id() << " at " << thePoint << endreq;  
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
    if ( m_pocaTool->ribbonDistances(*ribbon,track,*pocaData).isFailure() ) {
      MsgStream log(msgSvc(),name()) ;
      log << MSG::ERROR << "fillPocaData() call to ribbonDistances() failed" << endreq;
      return StatusCode::FAILURE;
    }
  } else if ( acdId.tile() ) {
    const AcdTileDim* tile = geomMap.getTile(acdId,*m_acdGeomSvc);
    if ( tile->statusCode().isFailure() ) {
      MsgStream log(msgSvc(),name()) ;
      log << MSG::ERROR << "Tile statusCode is failure" << endreq;
      return StatusCode::FAILURE;
    }
    if ( m_pocaTool->tileDistances(*tile,track,*pocaData).isFailure() ) {
      MsgStream log(msgSvc(),name()) ;
      log << MSG::ERROR << "fillPocaData() call to tileDistances() failed" << endreq;
      return StatusCode::FAILURE;
    }
  }  
  return StatusCode::SUCCESS;
}

StatusCode AcdTkrIntersectToolV2::holePoca(const AcdRecon::TrackData& /* track */, 
					   const AcdRecon::PocaData& /* pocaData */, const AcdTileDim& /* tile */,
					   std::vector<Event::AcdTkrGapPoca*>& /* gapPocas */) {

  return StatusCode::SUCCESS ;
}


StatusCode AcdTkrIntersectToolV2::gapPocaTile(const AcdRecon::TrackData& track, const AcdRecon::ExitData& /* data */,
					      const AcdRecon::PocaData& pocaData, AcdGeomMap& geomMap, 
					      std::vector<Event::AcdTkrGapPoca*>& gapPocas){
  
  signed char gapType = (unsigned char)(AcdRecon::None);
  unsigned char face = (unsigned char)(pocaData.m_id.face());
  unsigned char row = (unsigned char)(pocaData.m_id.row());
  unsigned char col = (unsigned char)(pocaData.m_id.column());
  unsigned char gap = 0;

  // bottow side tile have gaps only at the ends
  if ( face != 0 && row == 3 ) {
    gapType = AcdRecon::SideCornerEdge;
    if ( pocaData.m_inPlane.x() > 0 ) { gap = 1; }    
  } else {
    // if ( pocaData.m_inPlane.y() > 0 ) { row++; }    
    switch ( face ) {
    case 0: // TOP, Normal col and row counting 
      if ( pocaData.m_inPlane.x() > 0  && col == 4 ) {
	gapType = AcdRecon::TopCornerEdge;
      } else if ( pocaData.m_inPlane.x() < 0  && col == 0 ) {
	gapType = AcdRecon::TopCornerEdge;
      }	else {
	gapType = AcdRecon::Y_RibbonTop;    
      }
      break;
    case 1:  // -Y face, column go opposite of local x
    case 2:  // -X face, column go opposite of local x
      if ( pocaData.m_inPlane.x() < 0  && col == 4 ) {
	gapType = AcdRecon::SideCornerEdge;
      } else if ( pocaData.m_inPlane.x() > 0  && col == 0 ) {
	gapType = AcdRecon::SideCornerEdge; 
      } else {
	gapType = face == 1 ? AcdRecon::X_RibbonSide : AcdRecon::Y_RibbonSide;
      }
      break;
    case 3:  // +Y face, column go with local x
    case 4:
      if ( pocaData.m_inPlane.x() > 0 && col == 4 ) {
	gapType = AcdRecon::SideCornerEdge;
      } else if ( pocaData.m_inPlane.x() < 0 && col == 0 ) {
	gapType = AcdRecon::SideCornerEdge;
      } else {
	gapType = face == 3 ? AcdRecon::X_RibbonSide : AcdRecon::Y_RibbonSide;	
      }
      break;
    }
  }
  
  AcdTileDim* tile = const_cast<AcdTileDim*>(geomMap.getTile(pocaData.m_id,*m_acdGeomSvc));
  StatusCode sc = tile->latchGaps();
  if ( sc.isFailure() ) {
    return sc;
  }
  
  float gapSize = tile->gapAtLocalY(( pocaData.m_inPlane.x() < 0 ? -1 : 1), pocaData.m_inPlane.y());
  if ( gapSize < 0.001 ) {
    return sc;
  }

  idents::AcdGapId gapId(gapType,gap,face,row,col);

  Event::AcdTkrGapPoca* poca(0);
  sc = makeGapPoca(gapId,track,pocaData,gapSize,poca);
  if ( sc.isFailure() ) return sc;

  gapPocas.push_back(poca);

  return StatusCode::SUCCESS;
}


StatusCode AcdTkrIntersectToolV2::makeGapPoca(idents::AcdGapId& gapId, const AcdRecon::TrackData& track, 
					      const AcdRecon::PocaData& pocaData,
					      double gapSize, Event::AcdTkrGapPoca*& poca) {

  poca = 0;

  double arcLengthPlane = track.m_upward ? pocaData.m_arcLengthPlane : -1* pocaData.m_arcLengthPlane;
  double arcLengthPoca = track.m_upward ? pocaData.m_arcLength : -1* pocaData.m_arcLength;
  
  float vetoSigmaHit(0.);
  float local[2];
  local[0] = pocaData.m_activeX;

  float active[2];
  active[0] = pocaData.m_activeX; active[1] = pocaData.m_activeY;

  float halfGapSize = gapSize/2.;
  float gapCenter = -1.* halfGapSize;
  
  float gapNearDist(0.), gapFarDist(0.), centerGapDist(0.);
  
  if ( pocaData.m_activeX > gapCenter ) {
    // on tile side of gap center
    gapNearDist = pocaData.m_activeX;
    gapFarDist = gapNearDist + gapSize;
    centerGapDist = gapNearDist + halfGapSize;
  } else {
    // on far side of gap center
    gapNearDist = pocaData.m_activeX + gapSize;
    gapFarDist = pocaData.m_activeX;
    centerGapDist = gapNearDist + halfGapSize;   
  }
  
  // switch to active dist
  float gapActiveDist = halfGapSize - fabs(centerGapDist);

  switch ( gapId.gapType() ) {
  case AcdRecon::CornerRay:
  case AcdRecon::SideCornerEdge:
    // this gap has been calculated from a corner ray, no chance it could have a singal
    vetoSigmaHit = 0.;
    local[1] = pocaData.m_activeY;
    break;
  case AcdRecon::X_RibbonSide:
  case AcdRecon::Y_RibbonSide: 
  case AcdRecon::Y_RibbonTop: 
  default:
    // this gap has been calculated from a ribbon, this is some change in could have a signal
    vetoSigmaHit = 1.5;
    local[1] = ribbonLengthEstimate( gapId.face(), pocaData.m_hitsPlane );
    break;
   }

  // temp storage
  Point p( pocaData.m_poca.x(),pocaData.m_poca.y(), pocaData.m_poca.z());
  Vector voca( pocaData.m_voca.x(),pocaData.m_voca.y(), pocaData.m_voca.z());


  float activeDistErrProj = sqrt ( pocaData.m_planeError_proj(1,1) );
  float activeDistErrProp = sqrt ( pocaData.m_planeError_prop(1,1) );

  float gapNearSigmaProj = gapNearDist / activeDistErrProj;  
  float gapFarSigmaProj = gapFarDist / activeDistErrProj;

  float gapNearSigmaProp = gapNearDist / activeDistErrProp;  
  float gapFarSigmaProp = gapFarDist / activeDistErrProp;

  float vetoSigmaProj = sigmaEquivalent(gapNearSigmaProj,gapFarSigmaProj);
  float vetoSigmaProp = sigmaEquivalent(gapNearSigmaProp,gapFarSigmaProp);

  if ( false ) {
    std::cout << "Gap " << gapNearDist << ' ' << gapFarDist << ' '  << gapActiveDist << ' ' 
	      << gapNearSigmaProj << ' ' << gapFarSigmaProj << ' ' 
	      << vetoSigmaProj << std::endl 
	      << gapNearSigmaProp << ' ' << gapFarSigmaProp << ' ' 
	      << vetoSigmaProp << std::endl;
  }

  poca = new Event::AcdTkrGapPoca(gapId,track.m_index,active,				  
				  vetoSigmaHit,vetoSigmaProj,vetoSigmaProp,
				  pocaData.m_volumePlane,arcLengthPlane,pocaData.m_cosTheta,
				  pocaData.m_hitsPlane,local,pocaData.m_planeError_proj,pocaData.m_planeError_prop,
				  pocaData.m_volume,pocaData.m_region,arcLengthPoca,
				  gapActiveDist,activeDistErrProj,activeDistErrProp,
				  p,voca);
  
  return StatusCode::SUCCESS;
}

StatusCode AcdTkrIntersectToolV2::makeTkrPoint(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
					     Event::AcdTkrPoint*& tkrPoint ) {
  tkrPoint = 0;  
  float local[2];
  local[0] = data.m_x.x();
  local[1] = data.m_x.y();
  switch ( data.m_face ) {
  case 0:
  case 5:
    break;
  case 1:
  case 3:
    local[0] = data.m_x.y();
    local[1] = data.m_x.z();
    break;
  case 2:
  case 4:
    local[1] = data.m_x.z();
    break;    
  }
  tkrPoint = new Event::AcdTkrPoint(track.m_index,
				    data.m_face,data.m_arcLength,0.,
				    data.m_x,local,data.m_planeError_proj,data.m_planeError_prop);
  return StatusCode::SUCCESS;
}


float AcdTkrIntersectToolV2::sigmaEquivalent(float x1, float x2) {

  static const double root2 = sqrt(2.);
  double erf1 = TMath::Erf(x1);
  double erf2 = TMath::Erf(x2);
  double delErf = fabs(erf1 - erf2);
  if ( delErf < 1e-9 ) {
    // Far away, just return the nearer edge
    return fabs(x1 * root2);
  } else if ( delErf > 1. ) {
    // At least 50. prob to be from inside gap, 
    return 0.;
  }
 
  // Convert to the prob^2 that is inside the gap to sigma
  float sigmaSpan2 = -1. * log(delErf);

  // Check to see if we project inside the gap
  // (x1,x2) would be opposite signs
  if ( (x1 * x2) < 0 ) {
    return sqrt(sigmaSpan2);
  }

  // Not inside the gap, add sigmaSpan2 + x1*x1 (ie, in quadrature)
  sigmaSpan2 += (x1*x1);
  return sqrt(sigmaSpan2);
  
}

float AcdTkrIntersectToolV2::ribbonLengthEstimate(int face, const HepPoint3D& point) {

  static const AcdRecon::AcdVolume vol;

  float retVal(0.);

  switch ( face ) {
  case 0: // top: take the Y coord
    retVal = point.y();
    break;
  case 1: // -x side: take width of top + distance down side (negative values)
  case 2: // -y side: take width of top + distance down side (negative values)
    retVal -= vol.m_sides;
    retVal -= (vol.m_top - point.z());
    break;
  case 3: // +x side: take width of top + distance down side
  case 4: // +y side: take width of top + distance down side
    retVal += vol.m_sides;
    retVal += (vol.m_top - point.z());
  default:
    break;
  }

  return retVal;
}

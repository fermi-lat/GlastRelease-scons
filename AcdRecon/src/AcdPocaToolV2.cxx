#include "AcdPocaToolV2.h"

#include "../AcdRecon/AcdReconFuncsV2.h"

#include "AcdUtil/AcdTileDim.h"
#include "AcdUtil/AcdRibbonDim.h"
#include "AcdUtil/AcdUtilFuncs.h"

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
#include "geometry/Point.h"
#include "geometry/Vector.h"

#include "idents/AcdId.h"

DECLARE_TOOL_FACTORY(AcdPocaToolV2) ;

const double AcdPocaToolV2::MaxDoca(2000.);

AcdPocaToolV2::AcdPocaToolV2
 ( const std::string & type, 
   const std::string & name,
   const IInterface * parent )
 : AlgTool( type, name, parent )
 { 
   declareInterface<AcdIPocaToolV2>(this) ; 
   declareProperty("distanceCut",m_distanceCut=1999.);
   declareProperty("sigmaCut",m_sigmaCut=5.);
 }

AcdPocaToolV2::~AcdPocaToolV2()
{} 

// This function extracts geometry constants from xml
// file using GlastDetSvc
StatusCode AcdPocaToolV2::initialize()
 {
  MsgStream log(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;
  log<<MSG::INFO<<"BEGIN initialize()"<<endreq ;
 
  return sc;
 }



StatusCode AcdPocaToolV2::tileDistances (const AcdTileDim& tile,
					 const AcdRecon::TrackData& aTrack, 
					 AcdRecon::PocaData& data) {

  // initialize and sanity check
  data.reset(MaxDoca);
  StatusCode sc = tile.statusCode();

  MsgStream log(msgSvc(), name());
  
  if (sc.isFailure()) { 
    log << MSG::ERROR << "Tile Geometry data " << tile.acdId().id() << " is not ready." << endreq;
    return sc;
  }
  
  if ( log.level() <= MSG::DEBUG ) {
    log << MSG::DEBUG << "AcdHitPoca for track " << aTrack.m_index << ", tile " << tile.acdId().id() << endreq;  
  }

  // latch the tile id
  data.m_id = tile.acdId();
  
  // First, project track to plane of the tile
  HepPoint3D p;
  AcdRecon::ReconFunctions::tilePlane(aTrack,tile,data);

  // had to go backward to cross plane.  return
  if ( data.m_arcLengthPlane < 1e-9 ) return sc;

  bool isInside = data.m_active2D > 0;
  AcdRecon::ReconFunctions::tileEdgePoca(aTrack,tile,data,isInside);
  
  if ( log.level() <= MSG::DEBUG ) {    
    log << MSG::DEBUG << "AcdTilePlane.\n\t" 
	<< "s= " << data.m_arcLengthPlane 
	<< "; X_g=" << data.m_hitsPlane << "; X_l=" << data.m_inPlane << "; vol= " << data.m_volume << std::endl
	<< "\tcos(th)= " << data.m_cosTheta << "; l= " << data.m_path 
	<< "; active= (" << data.m_activeX << ',' << data.m_activeY << '|' << data.m_active2D << ')' << endreq;
    log << MSG::DEBUG << "AcdTilePoca.\n\t" 
	<< "s= " << data.m_arcLength << "; P= " << data.m_poca << "; V= " << data.m_voca 
	<< "; active3D= " << data.m_active3D << "; region= " << data.m_region << endreq;
  }
  return sc;
}


StatusCode AcdPocaToolV2::ribbonDistances(const AcdRibbonDim& ribbon,
					  const AcdRecon::TrackData& aTrack, 
					  AcdRecon::PocaData& data) {
  //Purpose and Method:  Calculate ActiveDistance for Ribbons
  // To simplify the situation - we assume the ribbons are merely lines 
  // Since the ribbons are segmented in the geometry - we treat the side segments and 
  // top segments separately, taking the minimum distance over all hit segmen

  // initialize and sanity check
  data.reset(MaxDoca);  // this is an active dist caculation, want to pick larger numbers
  StatusCode sc = ribbon.statusCode();

  MsgStream log(msgSvc(), name());  
  if ( log.level() <= MSG::DEBUG ) {
    log << MSG::DEBUG << "AcdHitPoca for track " << aTrack.m_index << ", ribbon  " << ribbon.acdId().id() << endreq;
  } 

  if (sc.isFailure()) {
    log << MSG::ERROR << "Ribbon Geometry data " << ribbon.acdId().id() << " is not ready." << endreq;
    return sc;
  }

  // latch the tile id
  data.m_id = ribbon.acdId();
  
  // Then project track to plane of the ribbon
  HepPoint3D p;

  // Then do the 3d poca to the ribbon
  AcdRecon::ReconFunctions::ribbonPoca(aTrack,ribbon,data);

  if ( log.level() <= MSG::DEBUG ) {    
    log << MSG::DEBUG << "AcdRibbonPlane.\n\t"
	<< "s= " << data.m_arcLengthPlane 
	<< "; X_g=" << data.m_hitsPlane << "; X_l=" << data.m_inPlane << "; vol= " << data.m_volume << std::endl
	<< "\tcos(th)= " << data.m_cosTheta << "; l= " << data.m_path 
	<< "; active= (" << data.m_activeX << ',' << data.m_activeY << '|' << data.m_active2D << ')' << endreq;  
    log << MSG::DEBUG << "AcdRibbonPoca.\n\t"
	<< "s= " << data.m_arcLength << "; R= " << data.m_rayLength 
	<< "; P=" << data.m_poca << "; V=" << data.m_voca 
	<< "; active= " << data.m_active3D << "; region= " << data.m_region << endreq;
  }
  return sc;
}


StatusCode AcdPocaToolV2::makePoca(const AcdRecon::TrackData& aTrack, 
				   const AcdRecon::PocaData& pocaData,	
				   const Event::AcdHit* theHit,			 
				   Event::AcdTkrHitPoca*& poca) {
  poca = 0;  
  const idents::AcdId& acdId = pocaData.m_id;

  MsgStream log(msgSvc(), name());  

  float local[2]; 
  float active2d[2];
  float mips[2];
  unsigned short flags[2];
  mips[0] =  theHit ? theHit->mips( Event::AcdHit::A) : 0.;
  mips[1] =  theHit ? theHit->mips( Event::AcdHit::B) : 0.;
  flags[0] =  theHit ? theHit->getFlags( Event::AcdHit::A) : 0;
  flags[1] =  theHit ? theHit->getFlags( Event::AcdHit::B) : 0;

  local[0] = pocaData.m_inPlane.x();
  active2d[0] = pocaData.m_activeX;
  active2d[1] = pocaData.m_activeY;
  float arcLengthPlane = aTrack.m_upward ? pocaData.m_arcLengthPlane : -1* pocaData.m_arcLengthPlane;
  float arcLengthPoca = aTrack.m_upward ? pocaData.m_arcLength : -1* pocaData.m_arcLength;
  float activeDist3DErrProj = pocaData.m_active3DErr_proj;
  float activeDist3DErrProp = pocaData.m_active3DErr_prop;

  float vetoSigmaProj(0.); 
  float vetoSigmaProp(0.); 

  if ( pocaData.m_active3DErr_prop < 1e-9 ) {
    activeDist3DErrProp = activeDist3DErrProj;
    if ( aTrack.m_index >= 0 ) {
      HepVector3D dx = aTrack.m_point - aTrack.m_current;    
      log << MSG::WARNING << "No Prop error for poca with " << acdId.id() << ' ' << aTrack.m_index << " at " 
          << pocaData.m_arcLengthPlane << ' ' << pocaData.m_arcLength << ' ' 
          << dx.mag() << ' ' << pocaData.m_active2D << ' ' <<  pocaData.m_active3D << std::endl;
    } else {
      // don't do sigma for vertices
      activeDist3DErrProj = -1.;
      activeDist3DErrProp = -1.;
    }    
  } 

  if ( activeDist3DErrProj > 0 ) {
    // We have a good projection error
    // Note that vetoSigmaProj can go negative if we are inside the tile or ribbon
    vetoSigmaProj = -1. * pocaData.m_active3D / activeDist3DErrProj;
  }

  if ( activeDist3DErrProp > 0 ) {
    // We have a good propagation error
    // Note that vetoSigmaProj can go negative if we are inside the tile or ribbon
    vetoSigmaProp = -1. * pocaData.m_active3D / activeDist3DErrProp;
  }

  float distance(0.);
  if ( acdId.tile() ) {
    // for tile local data is just active distances in X and Y, 
    // doca is active2D if in plane, otherwise active3D
    distance = pocaData.m_active2D > 0 ? pocaData.m_active2D : pocaData.m_active3D;
    local[1] = pocaData.m_inPlane.y();
  } else if ( acdId.ribbon() ) {    
    distance = pocaData.m_active3D;    
    local[1] = pocaData.m_ribbonLength;
  } else {
    // NA channel, do nothing
    return StatusCode::SUCCESS;
  }

  bool hasHit = (mips[0] > 0.001) || (mips[1] > 0.001);
  // default value in case we don't have a hit is 1.e4
  float vetoSigmaHit(1.e4);
  
  if ( hasHit ) {
    float totalMips = mips[0] + mips[1];
    // Note that vetoSigmaHit takes the difference between expected
    // and measured MIP signal and divides by the fit width of the
    // MIP peak at normal incidence (~0.45 for tiles; 1.0 for ribbons).
    // Path length corrections are made.
    if ( acdId.tile() ) {
      float expectedMips = 2. / pocaData.m_cosTheta;
      float sigmaMips = 0.45 / sqrt(fabs(pocaData.m_cosTheta));
      vetoSigmaHit = ( expectedMips - totalMips ) / sigmaMips;  
    } else {
      float expectedMips = 1.5;
      float sigmaMips = 1.;
      vetoSigmaHit = ( expectedMips - totalMips ) / sigmaMips;
    }
  }


  if ( false ) {
    log << MSG::WARNING << "Hit " << acdId.id() << ' ' << pocaData.m_active3D << ' ' 
        << activeDist3DErrProj << ' ' << activeDist3DErrProp << ' ' 
        << vetoSigmaHit << ' ' << vetoSigmaProj << ' ' << vetoSigmaProp << std::endl; 
  }  

  // temp storage
  Point p( pocaData.m_poca.x(),pocaData.m_poca.y(), pocaData.m_poca.z());
  Vector voca( pocaData.m_voca.x(),pocaData.m_voca.y(), pocaData.m_voca.z());

  poca = new Event::AcdTkrHitPoca(pocaData.m_id,aTrack.m_index,
				  active2d,mips,
				  vetoSigmaHit,vetoSigmaProj,vetoSigmaProp,
				  pocaData.m_volumePlane,arcLengthPlane,pocaData.m_cosTheta,
				  pocaData.m_hitsPlane,local,
				  pocaData.m_planeError_proj,pocaData.m_planeError_prop,
				  pocaData.m_volume,pocaData.m_region,arcLengthPoca,
				  distance,activeDist3DErrProj,activeDist3DErrProp,
				  p,voca,flags);

  return StatusCode::SUCCESS;

}


// @brief add a poca onto a list, subject to filtering cuts
StatusCode AcdPocaToolV2::filter(const AcdRecon::PocaDataMap& in, AcdRecon::PocaDataPtrMap& out) {
  StatusCode sc = StatusCode::SUCCESS;
  // loop on all the pocas we computed
  for ( AcdRecon::PocaDataMap::const_iterator it = in.begin(); it != in.end(); it++ ) {
    const AcdRecon::PocaData& data = it->second;
    if ( data.m_active3D < (-1.*m_distanceCut) ) {
      continue;
    }
    //if ( data.m_arcLengthPlane < -1e-5 ) continue;
    AcdRecon::PocaData* ptr = const_cast<AcdRecon::PocaData*>(&data);
    out[it->first] = ptr;
  }
  return sc;
}


// @@brief gets the amount of energy in cones 15,30 and 45 degrees around object direction
StatusCode AcdPocaToolV2::getConeDistances(const std::vector<Event::AcdTkrHitPoca*>& hitPocae, 
					   float& energy15, float& energy30, float& energy45,
					   float& triggerEnergy15, float& triggerEnergy30, float& triggerEnergy45) {
  bool ok = AcdUtil::UtilFunctions::getConeEnergies(hitPocae,
						    energy15,energy30,energy45,
						    triggerEnergy15,triggerEnergy30,triggerEnergy45);
  return StatusCode::SUCCESS;
}


// @brief filters out all but the best few POCA for a given object
StatusCode AcdPocaToolV2::selectPocae(std::vector<Event::AcdTkrHitPoca*>& hitPocae,
				      int nBest,
				      int nHitBest,
				      int nTrigBest) {

  std::list<Event::AcdTkrHitPoca*> toSave;
  std::list<Event::AcdTkrHitPoca*> toDelete;
  
  int nSave(0);
  int nHitSave(0);
  int nTrigSave(0);

  for ( std::vector<Event::AcdTkrHitPoca*>::iterator itr = hitPocae.begin(); itr != hitPocae.end(); itr++ ) {
    Event::AcdTkrHitPoca* hitPoca = *itr;
    bool save = false;
    if ( nBest == 0 || nSave < nBest ) save = true;
    if ( hitPoca->hasHit() && ( nHitBest == 0 || nHitSave < nHitBest ) ) save = true;
    if ( hitPoca->getTriggerVeto() && ( nTrigBest == 0 || nTrigSave < nTrigBest ) ) save = true;
    if ( save ) {
      nSave++;
      if ( hitPoca->hasHit() ) nHitSave++;
      if (  hitPoca->getTriggerVeto() ) nTrigSave++;
      toSave.push_back(hitPoca);
    } else {
      toDelete.push_back(hitPoca);
    }
  }
  // clear the vector
  hitPocae.clear();
  for ( std::list<Event::AcdTkrHitPoca*>::iterator itrS = toSave.begin(); itrS != toSave.end(); itrS++ ) {
    hitPocae.push_back(*itrS);
  }
  for ( std::list<Event::AcdTkrHitPoca*>::iterator itrD = toDelete.begin(); itrD != toDelete.end(); itrD++ ) {
    delete *itrD;
  }  
  return StatusCode::SUCCESS;  
}

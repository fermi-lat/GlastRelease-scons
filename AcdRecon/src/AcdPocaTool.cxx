#include "AcdPocaTool.h"

#include "../AcdRecon/AcdReconFuncs.h"

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
#include "../AcdRecon/RayDoca.h"
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
   declareProperty("distanceCut",m_distanceCut=1999.);
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



StatusCode AcdPocaTool::tileDistances (const AcdTileDim& tile,
				       const AcdRecon::TrackData& aTrack, 
				       AcdRecon::PocaData& data) {

  // initialize and sanity check
  data.reset(MaxDoca);
  StatusCode sc = tile.statusCode();

  MsgStream log(msgSvc(), name());
  
  if (sc.isFailure()) return sc;
  
  if ( log.level() <= MSG::DEBUG ) {
    log << MSG::DEBUG << "AcdHitPoca for track " << aTrack.m_index << ", tile " << tile.acdId().id() << endreq;  
  }

  // latch the tile id
  data.m_id = tile.acdId();
  
  // First, project track to plane of the tile
  HepPoint3D p;
  AcdRecon::tilePlane(aTrack,tile,data);

  // had to go backward to cross plane.  return
  if ( data.m_arcLengthPlane < 1e-9 ) return sc;

  if (data.m_active2D > 0) { 
    // get the distance to the closest edge
    AcdRecon::tileEdgePoca(aTrack,tile,data.m_arcLength,data.m_active3D,data.m_poca,data.m_pocaVector,data.m_region);
  } else {
    // get the distance to the closest edge or corner
    AcdRecon::tileEdgeCornerPoca(aTrack,tile,data.m_arcLength,data.m_active3D,data.m_poca,data.m_pocaVector,data.m_region);
  }
  
  if ( log.level() <= MSG::DEBUG ) {    
    log << MSG::DEBUG << "AcdTilePlane.\n\t" 
	<< "s= " << data.m_arcLengthPlane 
	<< "; X_g=" << data.m_hitsPlane << "; X_l=" << data.m_inPlane << "; vol= " << data.m_volume << std::endl
	<< "\tcos(th)= " << data.m_cosTheta << "; l= " << data.m_path 
	<< "; active= (" << data.m_activeX << ',' << data.m_activeY << '|' << data.m_active2D << ')' << endreq;
    log << MSG::DEBUG << "AcdTilePoca.\n\t" 
	<< "s= " << data.m_arcLength << "; P= " << data.m_poca << "; V= " << data.m_pocaVector 
	<< "; active3D= " << data.m_active3D << "; region= " << data.m_region << endreq;
  }
  return sc;
}


StatusCode AcdPocaTool::ribbonDistances(const AcdRibbonDim& ribbon,
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

  if (sc.isFailure()) return sc;

  // latch the tile id
  data.m_id = ribbon.acdId();
  
  // Then project track to plane of the ribbon
  HepPoint3D p;
  AcdRecon::ribbonPlane(aTrack,ribbon,data);

  // Then do the 3d poca to the ribbon
  AcdRecon::ribbonPoca(aTrack,ribbon,data.m_arcLength,data.m_ribbonLength,
		       data.m_active3D,data.m_poca,data.m_pocaVector,data.m_region);

  if ( log.level() <= MSG::DEBUG ) {    
    log << MSG::DEBUG << "AcdRibbonPlane.\n\t"
	<< "s= " << data.m_arcLengthPlane 
	<< "; X_g=" << data.m_hitsPlane << "; X_l=" << data.m_inPlane << "; vol= " << data.m_volume << std::endl
	<< "\tcos(th)= " << data.m_cosTheta << "; l= " << data.m_path 
	<< "; active= (" << data.m_activeX << ',' << data.m_activeY << '|' << data.m_active2D << ')' << endreq;  
    log << MSG::DEBUG << "AcdRibbonPoca.\n\t"
	<< "s= " << data.m_arcLength << "; R= " << data.m_ribbonLength 
	<< "; P=" << data.m_poca << "; V=" << data.m_pocaVector 
	<< "; active= " << data.m_active3D << "; region= " << data.m_region << endreq;
  }
  return sc;
}


StatusCode AcdPocaTool::makePoca(const AcdRecon::TrackData& aTrack, 
				 const AcdRecon::PocaData& pocaData,				 
				 Event::AcdTkrHitPoca*& poca) {
  poca = 0;  
  const idents::AcdId& acdId = pocaData.m_id;

  float local[2], distance(0.), arcLength(0.);
  if ( acdId.tile() ) {
    // for tile local data is just active distances in X and Y, 
    // doca is active2D if in plane, otherwise active3D
    local[0] = pocaData.m_activeX;
    local[1] = pocaData.m_activeY;
    distance = pocaData.m_active2D > 0 ? pocaData.m_active2D : pocaData.m_active3D;
    // This should probably be the same as the active distance (either 2D or 3D depending if it hit
    // the tile or not), except that there are pathological cases where the 3D cases is negative, because
    // it occurs behind the first point
    // Safer just to use the 2D one for now
    //   
    //double arcLength3D = aTrack.m_upward ? pocaData.m_arcLength : -1* pocaData.m_arcLength;
    double arcLengthPlane = aTrack.m_upward ? pocaData.m_arcLengthPlane : -1* pocaData.m_arcLengthPlane;
    arcLength = arcLengthPlane;
  } else if ( acdId.ribbon() ) {    
    // for ribbon local data is just active distances in local X (across ribbon) and 
    // total length along the ribbon
    // doca is active3D always
    local[0] = pocaData.m_active2D;
    local[1] = pocaData.m_ribbonLength;
    // This should probably be the same as the active distance (either 2D or 3D depending if it hit
    // the tile or not), except that there are pathological cases where the 3D cases is negative, because
    // it occurs behind the first point
    // Safer just to use the 2D one for now
    //   
    //double arcLength3D = aTrack.m_upward ? pocaData.m_arcLength : -1* pocaData.m_arcLength;
    double arcLengthPlane = aTrack.m_upward ? pocaData.m_arcLengthPlane : -1* pocaData.m_arcLengthPlane;
    arcLength = arcLengthPlane;
    distance = pocaData.m_active3D;    
  } else {
    // NA channel, do nothing
    return StatusCode::SUCCESS;
  }

  // temp storage
  static Event::AcdTkrLocalCoords localCoords;
  static Event::AcdPocaData pd;

  localCoords.set(local,pocaData.m_path,pocaData.m_cosTheta,pocaData.m_region,pocaData.m_planeError);
  pd.set(arcLength,distance,pocaData.m_active3DErr,pocaData.m_poca,pocaData.m_pocaVector);
  poca = new Event::AcdTkrHitPoca(pocaData.m_id,aTrack.m_index,localCoords,pd);

  return StatusCode::SUCCESS;

}


// @brief add a poca onto a list, subject to filtering cuts
StatusCode AcdPocaTool::filter(const AcdRecon::PocaDataMap& in, AcdRecon::PocaDataPtrMap& out) {
  StatusCode sc = StatusCode::SUCCESS;
  // loop on all the pocas we computed
  for ( AcdRecon::PocaDataMap::const_iterator it = in.begin(); it != in.end(); it++ ) {
    const AcdRecon::PocaData& data = it->second;
    if ( data.m_active3D < (-1.*m_distanceCut) ) continue;
    if ( data.m_arcLengthPlane < -1e-5 ) continue;

    AcdRecon::PocaData* ptr = const_cast<AcdRecon::PocaData*>(&data);
    out[it->first] = ptr;
  }
  return sc;
}

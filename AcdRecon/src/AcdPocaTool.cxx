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
   declareProperty("distanceCut",m_distanceCut=2000.);
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
  if (sc.isFailure()) return sc;

  // Get the Tile Center
  const HepPoint3D& xT = tile.tileCenter();  	
 
  // latch the tile id
  data.m_id = tile.acdId();

  // OK, first the POCA to the tile center
  HepPoint3D p;
  AcdRecon::pointPoca(aTrack,xT,data.m_arcLengthCenter,data.m_docaCenter,p);
  data.m_pocaCenter.set(p.x(),p.y(),p.z());

  // Then project track to plane of the tile
  AcdRecon::tilePlane(aTrack,tile,data.m_arcLengthPlane,data.m_localX,data.m_localY,
		      data.m_activeX,data.m_activeY,data.m_active2D,p);
  data.m_inPlane.set(p.x(),p.y(),p.z());

  if ( data.m_arcLengthPlane < 1e-9 ) return sc;

  if (data.m_active2D > 0) { 
    // get the distance to the closest edge
    AcdRecon::tileEdgePoca(aTrack,tile,data.m_arcLength,data.m_active3D,data.m_poca,data.m_pocaVector,data.m_region);
  } else {
    // get the distance to the closest edge or corner
    AcdRecon::tileEdgeCornerPoca(aTrack,tile,data.m_arcLength,data.m_active3D,data.m_poca,data.m_pocaVector,data.m_region);
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
  if (sc.isFailure()) return sc;

  // latch the tile id
  data.m_id = ribbon.acdId();
  
  // Then project track to plane of the ribbon
  HepPoint3D p;
  AcdRecon::ribbonPlane(aTrack,ribbon,data.m_arcLengthPlane,data.m_active2D,p);
  data.m_inPlane.set(p.x(),p.y(),p.z());

  // Then do the 3d poca to the ribbon
  AcdRecon::ribbonPoca(aTrack,ribbon,data.m_arcLength,data.m_active3D,data.m_poca,data.m_pocaVector,data.m_region);

  return sc;
}


StatusCode AcdPocaTool::makePoca(const AcdRecon::TrackData& aTrack, 
				 const AcdRecon::PocaData& pocaData,				 
				 Event::AcdTkrHitPoca*& poca) {
  poca = 0;
  
  double arcLength = aTrack.m_upward ? pocaData.m_arcLength : -1* pocaData.m_arcLength;

  idents::AcdId acdId = pocaData.m_id;
    
  float local[2];
  local[0] = pocaData.m_activeX;
  local[1] = pocaData.m_activeY;
  HepMatrix localCov(2,2);
  localCov[0][0] = pocaData.m_localCovXX;
  localCov[1][1] = pocaData.m_localCovYY;
  localCov[0][1] = localCov[1][0] = pocaData.m_localCovXY;
  float distance = pocaData.m_active2D > 0 ? pocaData.m_active2D : pocaData.m_active3D;

  // temp storage
  static Event::AcdTkrLocalCoords localCoords;
  static Event::AcdPocaData pd;

  localCoords.set(local,pocaData.m_path,pocaData.m_cosTheta,pocaData.m_region,localCov);
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
    if ( data.m_arcLength < 1e-5 ) continue;
    AcdRecon::PocaData* ptr = const_cast<AcdRecon::PocaData*>(&data);
    out[it->first] = ptr;
  }
  return sc;
}

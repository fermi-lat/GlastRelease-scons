
#ifndef __AcdITkrIntersectTool_H
#define __AcdITkrIntersectTool_H 1

#include "GaudiKernel/IAlgTool.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrTrackParams.h"
#include "Event/Recon/AcdRecon/AcdTkrIntersection.h"
#include "Event/Recon/AcdRecon/AcdTkrGapPoca.h"
#include "Event/Recon/AcdRecon/AcdTkrPoint.h"
#include <vector>

#include "../AcdRecon/AcdReconStruct.h"

class IPropagator;
class AcdGeomMap;

/**   
 * @class AcdITkrIntersectTool
 *
 * @brief Gaudi interface for Tool that uses Kalman propagator to extend tracks 
 * through GEANT model of detector
 *
 * Much of the actual calculations live in AcdRecon::AcdReconFuncs
 * This code primarily loops over the objects and calls the relavent calculations
 *
 * @author Eric Charles
 *
 * $Header$
 **/

static const InterfaceID IID_AcdITkrIntersectTool("AcdITkrIntersectTool",1,0) ;

class AcdITkrIntersectTool : virtual public IAlgTool {


 public:
  
  /// retrieve Gaudi interface ID
  static const InterfaceID& interfaceID()
    { return IID_AcdITkrIntersectTool; }
  
  AcdITkrIntersectTool() {}
  virtual ~AcdITkrIntersectTool() {}
  
  /** @brief propagate a track through the GEANT model
   *
   *  @param prop the kalman propagator (contains track parameters)
   *  @param track a few ancillary parameters (like direction, starting point)
   *  @param data where the track exists the ACD volume ( a cube around the ACD )
   *  @param pocaMap map containing the AcdRecon::PocaData objects by AcdId
   *  @param hitMap map containing mask of tile discriminator values
   *  @param geomMap map containing geometrical discriptions of tiles and ribbons by AcdId
   *  @param intersections TDS object contains the intersection we make
   *  @param gapPocas TDS object contains data about how close we came to gaps & ribbons
   *  @return Success or Failure
   **/
  virtual StatusCode makeIntersections(IPropagator& prop,
				       const AcdRecon::TrackData& track,
				       const AcdRecon::ExitData& data,	
				       const AcdRecon::PocaDataPtrMap& pocaMap,
				       const AcdRecon::AcdHitMap& hitMap,
				       AcdGeomMap& geomMap,
				       Event::AcdTkrIntersectionCol& intersections,
				       Event::AcdTkrGapPocaCol& gapPocas) = 0;
  
  /** @brief make the TDS object that states where the track left the ACD
   *
   *  @param track a few ancillary parameters (like direction, starting point)
   *  @param data where the track exists the ACD volume ( a cube around the ACD )
   *  @param params track parameters at exit point
   *  @param tkrPoint TDS object we make
   *  @return Success or Failure
   **/
  virtual StatusCode makeTkrPoint(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
				  const Event::TkrTrackParams& params, Event::AcdTkrPoint*& tkrPoint ) = 0;


} ;

#endif




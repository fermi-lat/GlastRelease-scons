
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
* Base class for clustering tools
*
* $Header$
*/

static const InterfaceID IID_AcdITkrIntersectTool("AcdITkrIntersectTool",1,0) ;

class AcdITkrIntersectTool : virtual public IAlgTool {


 public:
  
  // retrieve Gaudi interface ID
  static const InterfaceID& interfaceID()
    { return IID_AcdITkrIntersectTool; }
  
  AcdITkrIntersectTool() {}
  virtual ~AcdITkrIntersectTool() {}
  
  // @brief main method
  virtual StatusCode makeIntersections(IPropagator& prop,
				       const AcdRecon::TrackData& track,
				       const AcdRecon::ExitData& data,	
				       const AcdRecon::PocaDataPtrMap& pocaMap,
				       const AcdRecon::AcdHitMap& hitMap,
				       AcdGeomMap& geomMap,
				       Event::AcdTkrIntersectionCol& intersections,
				       Event::AcdTkrGapPocaCol& gapPocas) = 0;

  // @brief calculate the arclength at which a track exits the tracking volume
  virtual StatusCode exitsLAT(const Event::TkrTrack& track, bool forward,
			      AcdRecon::ExitData& data) = 0;

  // @brief make the TDS object that states where the track left the ACD
  virtual StatusCode makeTkrPoint(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
				  const Event::TkrTrackParams& params, Event::AcdTkrPoint*& tkrPoint ) = 0;


} ;

#endif




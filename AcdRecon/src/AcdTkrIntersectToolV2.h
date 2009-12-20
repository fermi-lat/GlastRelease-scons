
#ifndef __AcdTkrIntersectToolV2_H
#define __AcdTkrIntersectToolV2_H 1

#include "AcdITkrIntersectToolV2.h"
#include "AcdIPocaToolV2.h"
#include "AcdUtil/AcdGeomMap.h"

#include "../AcdRecon/AcdReconStruct.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "AcdUtil/IAcdGeometrySvc.h"
#include "GlastSvc/Reco/IPropagator.h" 
#include "TkrUtil/ITkrGeometrySvc.h"

#include "idents/VolumeIdentifier.h"

#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrTrackParams.h"



class MsgStream;
namespace CLHEP {class HepMatrix;}

/**   
 * @class AcdTkrIntersectToolV2
 *
 * @brief Gaudi Tool that uses Kalman propagator to extend tracks through GEANT model 
 * of detector
 *
 * Much of the actual calculations live in AcdRecon::AcdReconFuncs
 * This code primarily loops over the objects and calls the relavent calculations
 *
 * Inputs:
 *  - Called for every track in both directions
 *    - Kalman propagator
 *    - Track direction, starting point
 *    - Data about where track leaves ACD volume
 *    - Map of all existing POCA calculations for this track, by AcdId
 *    - Map of discrimator status for all ACD channels
 *    - Acd geometry map
 *
 * TDS Outputs:
 *  - Event::AcdTkrIntersectionCol: all the intersections with GEANT model 
 *  - Event::AcdTkrGapPocaCol: POCA w.r.t. gaps and ribbons in the GEANT model 
 *  - Event::AcdTkrPointCol: Data about where the track exits the ACD volume
 *
 * Algorithm:
 *  - Phase I: Called for every track in both directions,
 *    Called with propagator, track exit point, map of already calculated POCAs, map of hit data 
 *    propagator has already been fed final arclength during initialization
 *    - Sort POCAs by arclength
 *    - run track out to final arclength
 *    - loop on propagator steps, ignoring all non-ACD elements
 *      - check to see if any POCAs have been crossed
 *      - fill AcdRecon::PocaData for any new elements (ie, w/ no signal)
 *      - use track parameters to calculate error for any POCAs we crossed
 *      - latch tile and ribbons crossings for gap calculations later
 *      - build Event::AcdTkrIntersection object
 *   - build Event::AcdTkrGapPoca object storing data to nearest gap
 *     - if we hit a ribbon use that with gapPocaRibbon()
 *     - else if we hit a tile use that with gapPocaTile()
 *     - else if we didn't go out bottom of ACD use fallbackToNominal()
 *   - build Event::AcdTkrGapPoca to corner edge rays

 *
 * This tool has no JO parameters
 +
 * @author Eric Charles
 *
 * $Header$
*/


class AcdTkrIntersectToolV2 : public AcdITkrIntersectToolV2,  public AlgTool {
	
 public:
    
  AcdTkrIntersectToolV2
    ( const std::string & type, 
      const std::string & name,
      const IInterface * parent ) ;
  virtual ~AcdTkrIntersectToolV2() ;
  
  /// @brief Intialization of the tool
  virtual StatusCode initialize();

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
				       const Event::TkrTrackParams& trackParams,
				       const AcdRecon::TrackData& track,
				       const AcdRecon::ExitData& data,	
				       const AcdRecon::AcdHitMap& hitMap,
				       AcdGeomMap& geomMap,
				       AcdRecon::PocaDataPtrMap& pocaMap,
				       std::list<AcdRecon::PocaData> ownedPocaData,
				       std::vector<Event::AcdTkrGapPoca*>& gapPocas);


  /** @brief make the TDS object that states where the track left the ACD
   *
   *  @param track a few ancillary parameters (like direction, starting point)
   *  @param data where the track exists the ACD volume ( a cube around the ACD )
   *  @param tkrPoint TDS object we make
   *  @return Success or Failure
   **/
  virtual StatusCode makeTkrPoint(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
				  Event::AcdTkrPoint*& tkrPoint );


protected:
  

  
  /** @brief fill extra data about poca when once we have propagated track to an element
   *
   *  This is called when the tile or ribbon doesn't have a signal, but was crossed 
   *  during the GEANT propagation phase
   *
   *  @param track projection data
   *  @param volId GEANT volume in question
   *  @param pocaMap map containing the AcdRecon::PocaData objects by AcdId
   *  @param arc arclength at which the poca occurs
   *  @param geomMap map containing geometrical discriptions of tiles and ribbons by AcdId
   *  @param data where the track exists the ACD volume ( a cube around the ACD )
   *  @param ownedPocaData poca data objects made by this tool
   *  @param pocaData object with POCA calculation data
   *  @return Success or Failure
   **/
  virtual StatusCode fillPocaData(const AcdRecon::TrackData& track, const idents::VolumeIdentifier& volId,
				  const AcdRecon::PocaDataPtrMap& pocaMap, const double& arc,
				  AcdGeomMap& geomMap,
				  std::list<AcdRecon::PocaData>& ownedPocaData, AcdRecon::PocaData*& pocaData);

  /** @brief fill data about how close we came to a tile screw hole
   *
   *  @param track projection data
   *  @param pocaData object with POCA calculation data
   *  @param tile description of tile geometry
   *  @param gapPocas TDS object contains data about how close we came to gaps (tile hole in this case)
   *  @return Success or Failure
   **/  
  virtual StatusCode holePoca(const AcdRecon::TrackData& track, const AcdRecon::PocaData& pocaData, const AcdTileDim& tile,
			      std::vector<Event::AcdTkrGapPoca*>& gapPocas);  
    
  /** @brief fill data about how close we came to a gap (could be ribbon or edge)
   *
   *  This is called when we crossed a tile during the GEANT propagation phase
   *
   *  @param track projection data
   *  @param data where the track exists the ACD volume ( a cube around the ACD )
   *  @param pocaData object with POCA calculation data
   *  @param geomMap map containing geometrical discriptions of tiles and ribbons by AcdId
   *  @param gapPocas TDS object contains data about how close we came to gaps (tile hole in this case)
   *  @return Success or Failure
   **/  
  virtual StatusCode gapPocaTile(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
				 const AcdRecon::PocaData& pocaData, AcdGeomMap& geomMap, 
				 std::vector<Event::AcdTkrGapPoca*>& gapPocas);


  /** @brief make a TDS object to describe how close we came to a gap
   *
   *  @param gapId, which gap
   *  @param track projection data
   *  @param pocaData about what we hit near the gap
   *  @param distance to gap   
   *  @param poca TDS object contains data about how close we came to gaps
   *  @return Success or Failure
   **/  
  virtual StatusCode makeGapPoca(idents::AcdGapId& gapId, const AcdRecon::TrackData& track, const AcdRecon::PocaData& pocaData,
				 double distance, Event::AcdTkrGapPoca*& poca);

  /// Check if a volume id is active ACD dectector
  static bool checkVolId(idents::VolumeIdentifier& volId);

  /// Get the mask for a given hit
  static unsigned char getHitMask(const Event::AcdHit& aHit);

  /// Get the sigma equivalent of an interval
  static float sigmaEquivalent(float x1, float x2);

  /// Estimate ribbon length
  static float ribbonLengthEstimate(int face, const HepPoint3D& point);

private:

  /// Tool to do poca calculations
  AcdIPocaToolV2*  m_pocaTool;
  /// G4 Propagator
  IPropagator *    m_G4PropTool; 
  /// Detector service, for handling volume IDs.
  IGlastDetSvc*    m_detSvc; 
  /// ACD geometry service
  IAcdGeometrySvc* m_acdGeomSvc;
 
} ;

#endif
	
	
	

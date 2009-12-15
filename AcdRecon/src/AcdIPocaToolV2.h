
#ifndef __AcdIPocaToolV2_H
#define __AcdIPocaToolV2_H 1

#include "GaudiKernel/IAlgTool.h"

#include "Event/Recon/AcdRecon/AcdHit.h"
#include "Event/Recon/AcdRecon/AcdTkrHitPoca.h"
#include "Event/Recon/AcdRecon/AcdPocaMap.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"

#include "../AcdRecon/AcdReconStruct.h"

namespace Event {
  class AcdHit;
}
class AcdTileDim;
class AcdRibbonDim;

/**   
 * @class AcdIPocaToolV2
 *
 * @brief Gaudi interface for Tool that calculates Point of Closest Approach (POCA) between 
 * ACD elements and track projections.
 *
 * The actual calculations live in AcdRecon::AcdReconFuncs
 * This code just loops over the objects and calls the relavent calculations
 *
 * @author Eric Charles
 *
 * $Header$
*/

static const InterfaceID IID_AcdIPocaToolV2("AcdIPocaToolV2",1,1) ;

class AcdIPocaToolV2 : virtual public IAlgTool {

public:
  
  /// retrieve Gaudi interface ID
  static const InterfaceID& interfaceID()
    { return IID_AcdIPocaToolV2; }
  
  AcdIPocaToolV2() {}
  virtual ~AcdIPocaToolV2() {}
    
  /** @brief calculate the distance of closest approach between the track and the tile
   *
   *   This includes the distance of closest approach to the center of the tile
   *   and both the 2d and 3d distances to the closest edge or corner
   *
   *  @param tile geometrical description of the tile
   *  @param aTrack track projection data
   *  @param data object to cache output of calculations
   *  @return Success or Failure
   **/
  virtual StatusCode tileDistances (const AcdTileDim& tile,
				    const AcdRecon::TrackData& aTrack, 
				    AcdRecon::PocaData& data) = 0;
  
  /** @brief calculate the distance of closest approach between the track and the ribbon
   *
   *   This includes both the 2d and 3d distances to the ray defining the ribbon segement
   *
   *  @param ribbon geometrical description of the ribbon
   *  @param aTrack track projection data
   *  @param data object to cache output of calculations
   *  @return Success or Failure
   **/
  virtual StatusCode ribbonDistances(const AcdRibbonDim& ribbon,
				     const AcdRecon::TrackData& aTrack, 
				     AcdRecon::PocaData& data) = 0;

  /** @brief Make an Event::AcdTrkHitPoca object, given the PocaData
   *
   *  @param aTrack track projection data
   *  @param data cached output of calculations
   *  @param poca newly made TDS object
   *  @return Success or Failure
   **/
  virtual StatusCode makePoca(const AcdRecon::TrackData& aTrack, 
			      const AcdRecon::PocaData& data,			      
			      const Event::AcdHit* theHit,
			      Event::AcdTkrHitPoca*& poca) = 0;

   /** @brief put the pocas into the output map, subject to filtering cuts
   *
   *  @param in all the caculated POCA
   *  @param only selected POCA objects
   *  @return Success or Failure
   **/ 
  virtual StatusCode filter(const AcdRecon::PocaDataMap& in, AcdRecon::PocaDataPtrMap& out) = 0;

} ;

#endif





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
   *  @param out only selected POCA objects
   *  @return Success or Failure
   **/ 
  virtual StatusCode filter(const AcdRecon::PocaDataMap& in, AcdRecon::PocaDataPtrMap& out) = 0;


   /** @brief gets the amount of energy in cones 15,30 and 45 degrees around object direction
   *
   *  @param hitPocae all the caculated POCA
   *  @param energy15 -> Energy within 15 degrees of the object
   *  @param energy30 -> Energy within 30 degrees of the object
   *  @param energy45 -> Energy within 45 degrees of the object
   *  @param triggerEnergy15 -> Energy within 15 degrees of the object from triggered tiles
   *  @param triggerEnergy30 -> Energy within 30 degrees of the object from triggered tiles
   *  @param triggerEnergy45 -> Energy within 45 degrees of the object from triggered tiles
   *  @return Success or Failure
   **/ 
  virtual StatusCode getConeDistances(const std::vector<Event::AcdTkrHitPoca*>& hitPocae, 
				      float& energy15, float& energy30, float& energy45,
				      float& triggerEnergy15, float& triggerEnergy30, float& triggerEnergy45) = 0;

  
  /** @brief filters out all but the best few POCA for a given object
   *
   *  @param hitPocae all the caculated POCA.  Note that this function will remove and delete the POCAe that are not selected
   *  @param nBest -> The number of poca to save (0 means all)
   *  @param nHitBest -> The number of poca for hit tiles or ribbon to save (0 means all)
   *  @param nTrigBest -> The number of poca for triggered tiles or ribbon to save (0 means all)
   *  @return Success or Failure
   **/ 
  virtual StatusCode selectPocae(std::vector<Event::AcdTkrHitPoca*>& hitPocae,
				 int nBest,
				 int nHitBest,
				 int nTrigBest) = 0;
				 
 
  

} ;

#endif




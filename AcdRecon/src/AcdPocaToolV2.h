
#ifndef __AcdPocaToolV2_H
#define __AcdPocaToolV2_H 1

#include "AcdIPocaToolV2.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ISvcLocator.h"

#include "CLHEP/Geometry/Point3D.h"
#include <vector>

class AcdTileDim;
class AcdRibbonDim;
class MsgStream;
class IGlastDetSvc;

class Point;
class Vector;

/**   
 * @class AcdPocaToolV2
 *
 * @brief Gaudi tool that calculates Point of Closest Approach (POCA) between 
 * ACD elements and track projections.
 *
 * The actual calculations live in AcdRecon::AcdReconFuncs
 * This code just loops over the objects and calls the relavent calculations
 *
 * Inputs:
 *  - For every track X tile/ribbon combination
 *    - AcdRecon::TrackData: track projection data
 *    - AcdTileDim or AcdRibbonDim: detector geometry data
 *
 * Intermediate Outputs:
 *  - AcdRecon::PocaDataMap: Poca Data from every track X tile/ribbon combination
 *  - AcdRecon::PocaDataPtrMap: Poca Data subjected to cuts 
 *
 * TDS Outputs:
 *  - Event::AcdTkrHitPoca object for every accepted track X tile/ribbon combination 
 *    - Stored in Event/AcdRecon/AcdTkrHitPocaCol collection
 *
 * Algorithm:
 *  - Phase I: Called for every track X tile/ribbon combination to calculate poca
 *    - For tiles call is to tileDistances()
 *      - Project track to tile plane using AcdRecon::tilePlane().  
 *        If intersection is in backward direction, stop and return.
 *        Otherwise this fills:
 *        - AcdRecon::PocaData::m_activeX
 *        - AcdRecon::PocaData::m_activeY
 *        - AcdRecon::PocaData::m_active2D
 *        - AcdRecon::PocaData::m_arcLengthPlane
 *        - AcdRecon::PocaData::m_hitsPlane
 *        - AcdRecon::PocaData::m_volume
 *      - Get 3D Poca to tile edge. This fills:
 *        - AcdRecon::PocaData::m_arcLength
 *        - AcdRecon::PocaData::m_active3D
 *        - AcdRecon::PocaData::m_poca
 *        - AcdRecon::PocaData::m_pocaVector
 *        - AcdRecon::PocaData::m_region
 *          - This uses AcdRecon::tileEdgePoca() if inside tile (m_active2D > 0) 
 *          - Or AcdRecon::tileEdgeCornerPoca() if outside tile
 *      
 *    - For ribbons call is to ribbonDistances():
 *      - Project track to ribbon plane (side or top of ACD) using AcdRecon::ribbonPlane().
 *        If intersection is in backward direction, stop and return.
 *        Otherwise this fills:
 *        - AcdRecon::PocaData::m_active2D
 *        - AcdRecon::PocaData::m_arcLengthPlane
 *        - AcdRecon::PocaData::m_hitsPlane
 *        - AcdRecon::PocaData::m_volume
 *      - Get 3D Poca to ribbon. This fills:
 *        - AcdRecon::PocaData::m_arcLength
 *        - AcdRecon::PocaData::m_ribbonLength
 *        - AcdRecon::PocaData::m_active3D
 *        - AcdRecon::PocaData::m_poca
 *        - AcdRecon::PocaData::m_pocaVector
 *        - AcdRecon::PocaData::m_region
 *
 *  - Phase II, Handed AcdRecon::PocaDataMap containing every AcdRecon::PocaData object.
 *    - Filters out all POCA with m_active3D < m_distanceCut.
 *
 *  - Phase III, Called for remaining AcdRecon::PocaData objects to make Event::AcdTkrHitPoca TDS objects
 *    - Set arcLength negative for downgoing intersections
 *    - For tiles:
 *      - Use PocaData::m_activeX, PocaData::m_activeY for local frame distances
 *      - For doca use PocaData::m_active2D if inside tile, PocaData::m_active3D otherwise
 *    - For ribbons:
 *      - Use PocaData::m_active2D, PocaData::m_ribbonLength for local frame distances
 *      - For doca use PocaData::m_active3D
 *    
 *
 * This tool has 2 JO paramters:  
 *  - distanceCut [1999.]   : Filter out POCA when doca is > this value
 *  - sigmaCut [5.]         : Unused!!  (Filter out POCA when doca/docaError) is > this value
 *
 *
 * @author Eric Charles
 *
 * $Header$
 **/


class AcdPocaToolV2 : public AcdIPocaToolV2,  public AlgTool {

public:
  
  static const double MaxDoca;
  
public:
    
  // @brief Standard Gaudi constructor, defines and sets parameters to default values
  AcdPocaToolV2
  ( const std::string & type, 
    const std::string & name,
    const IInterface * parent ) ;

  // @brief Trivial destructor
  virtual ~AcdPocaToolV2() ;
  
  // @brief Intialization of the tool
  virtual StatusCode initialize();
  
  // @brief calculate the distance of closest approach between the track and the tile
  //   This includes the distance of closest approach to the center of the tile
  //   and both the 2d and 3d distances to the closest edge or corner
  virtual StatusCode  tileDistances(const AcdTileDim& tile,
				    const AcdRecon::TrackData& aTrack, 
				    AcdRecon::PocaData& data);
  
  // @brief calculate the distance of closest approach between the track and the ribbon
  //   This includes both the 2d and 3d distances to the ray defining the ribbon segement
  virtual StatusCode ribbonDistances(const AcdRibbonDim& ribbon,
				     const AcdRecon::TrackData& aTrack, 
				     AcdRecon::PocaData& data);

  // @brief Make an AcdTrkPoca object, given the PocaData and the G4Propagator
  virtual StatusCode makePoca(const AcdRecon::TrackData& aTrack, 
			      const AcdRecon::PocaData& data, 
			      const Event::AcdHit* theHit,
			      Event::AcdTkrHitPoca*& poca);

  // @brief put the pocas onto a list, subject to filtering cuts
  virtual StatusCode filter(const AcdRecon::PocaDataMap& in, AcdRecon::PocaDataPtrMap& out);
  

  /** @brief gets the amount of energy in cones 15,30 and 45 degrees around object direction
   *
   *  @hitPocae all the caculated POCA
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
				      float& triggerEnergy15, float& triggerEnergy30, float& triggerEnergy45);


  
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
				 int nTrigBest);
				 

private:

  // parameters
  float m_distanceCut;
  float m_sigmaCut;

} ;

#endif
	
	
	


#ifndef __AcdIPocaTool_H
#define __AcdIPocaTool_H 1

#include "GaudiKernel/IAlgTool.h"

#include "Event/Recon/AcdRecon/AcdHit.h"
#include "Event/Recon/AcdRecon/AcdTkrHitPoca.h"
#include "Event/Recon/AcdRecon/AcdPocaMap.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"

#include "../AcdRecon/AcdReconStruct.h"

class AcdTileDim;
class AcdRibbonDim;

/**   
* @class AcdIPocaTool
*
* Base class for clustering tools
*
* $Header$
*/

static const InterfaceID IID_AcdIPocaTool("AcdIPocaTool",1,0) ;

class AcdIPocaTool : virtual public IAlgTool {

public:
  
  // retrieve Gaudi interface ID
  static const InterfaceID& interfaceID()
    { return IID_AcdIPocaTool; }
  
  AcdIPocaTool() {}
  virtual ~AcdIPocaTool() {}
    
  // @brief calculate the distance of closest approach between the track and the tile
  //   This includes the distance of closest approach to the center of the tile
  //   and both the 2d and 3d distances to the closest edge or corner
  virtual StatusCode tileDistances (const AcdTileDim& tile,
				    const AcdRecon::TrackData& aTrack, 
				    AcdRecon::PocaData& data) = 0;
  
  // @brief calculate the distance of closest approach between the track and the ribbon
  //   This includes both the 2d and 3d distances to the ray defining the ribbon segement
  virtual StatusCode ribbonDistances(const AcdRibbonDim& ribbon,
				     const AcdRecon::TrackData& aTrack, 
				     AcdRecon::PocaData& data) = 0;

  // @brief Make an AcdTrkPoca object, given the PocaData 
  virtual StatusCode makePoca(const AcdRecon::TrackData& aTrack, 
			      const AcdRecon::PocaData& data,			      
			      Event::AcdTkrHitPoca*& poca) = 0;

  // @brief put the pocas onto a list, subject to filtering cuts
  virtual StatusCode filter(const AcdRecon::PocaDataMap& in, AcdRecon::PocaDataPtrMap& out) = 0;

} ;

#endif




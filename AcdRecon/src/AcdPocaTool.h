
#ifndef __AcdPocaTool_H
#define __AcdPocaTool_H 1

#include "AcdIPocaTool.h"

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
* @class AcdPocaTool
*
* Class for calculation the distances between tracks and hit acd tiles and ribbons
*
* $ $
*/


class AcdPocaTool : public AcdIPocaTool,  public AlgTool {

public:
  
  static const double MaxDoca;
  
public:
    
  // @brief Standard Gaudi constructor, defines and sets parameters to default values
  AcdPocaTool
  ( const std::string & type, 
    const std::string & name,
    const IInterface * parent ) ;

  // @brief Trivial destructor
  virtual ~AcdPocaTool() ;
  
  // @brief Intialization of the tool
  virtual StatusCode initialize();
  
  // @brief calculate the distance of closest approach between the track and the tile
  //   This includes the distance of closest approach to the center of the tile
  //   and both the 2d and 3d distances to the closest edge or corner
  virtual StatusCode tileDistances (const AcdTileDim& tile,
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
			      Event::AcdTkrHitPoca*& poca);

  // @brief put the pocas onto a list, subject to filtering cuts
  virtual StatusCode filter(const AcdRecon::PocaDataMap& in, AcdRecon::PocaDataPtrMap& out);
  
private:

  // parameters
  float m_distanceCut;
  float m_sigmaCut;

} ;

#endif
	
	
	

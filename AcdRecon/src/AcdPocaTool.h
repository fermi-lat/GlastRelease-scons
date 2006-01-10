
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
class IPropagator;

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
  
  // @brief calculate the distance of closest approach between the track and the tile center
  virtual StatusCode doca (const AcdTileDim& tile,
			   const Event::TkrTrack& aTrack, 
			   PocaData& data);
  
  // @brief calculate the distance in the plane of the tile to the nearest edge of the tile
  //  this calculation is positive if the track hits the tile and negative otherwise
  //  the returned value is the largest number for any of the four tile edges
  virtual StatusCode hitTileDist(const AcdTileDim& tile,
				 const Event::TkrTrack& aTrack, 
				 PocaData& data);

  // @brief calculate the 3D distance to the nearest edge of the tile.  This calculation is 
  //  positive if the track hits the tile and negative otherwise
  //  the returned value is the largest number for any of the tile edges or corners
  virtual StatusCode tileActiveDist(const AcdTileDim& tile,
				    const Event::TkrTrack& aTrack, 
				    PocaData& data);
  
  // @brief calculate the 3D distance ribbon.  This calculation is 
  //  positive if the track hits the ribbon and negative otherwise
  virtual StatusCode hitRibbonDist(const AcdRibbonDim& ribbon,
				   const Event::TkrTrack& aTrack, 
				   PocaData& data);

  // @brief Make an AcdTrkPoca object, given the PocaData and the G4Propagator
  virtual StatusCode makePoca(const Event::TkrTrack& track, int iTrack,
			      const PocaData& poca, const idents::AcdId& acdId,
			      IPropagator& g4PropTool, Event::AcdTkrPoca*& pocaCol);

protected:

  // @brief calculate the 3D distance to the nearest edge or corner of the tile when the track 
  //  does not hit the tile.  This value should always be negative
  virtual StatusCode docaActiveDist(const AcdTileDim& tile,
				    const Event::TkrTrack& aTrack, 
				    PocaData& data);

  // @brief calculate the track parameters at the POCA
  virtual StatusCode getParamsAtPoca(const Event::TkrTrack& aTrack, bool forward, double arcLength,
				     IPropagator& g4PropTool,
				     Event::TkrTrackParams& paramsAtPoca);
  
  // @brief project the track error along the doca line
  virtual StatusCode projectError(const PocaData& poca, const Event::TkrTrackParams& paramsAtPoca,
				  double& pocaError);  

private:

  float m_distanceCut;
  float m_sigmaCut;

} ;

#endif
	
	
	

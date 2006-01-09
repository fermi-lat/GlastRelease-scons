
#ifndef __AcdIPocaTool_H
#define __AcdIPocaTool_H 1

#include "GaudiKernel/IAlgTool.h"

#include "Event/Recon/AcdRecon/AcdHit.h"
#include "Event/Recon/AcdRecon/AcdTkrPoca.h"
#include "Event/Recon/AcdRecon/AcdPocaMap.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"


class IPropagator;
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

  // This struct stores the data about the point of closest approach that we pass around internally
  struct PocaData {
  public:
    PocaData()
      :m_dist(0.),m_arcLength(-1.),m_region(0){;}
    void reset(double maxDist) {
      m_dist = maxDist;
      m_arcLength = -1;
      m_region = 0;
      m_poca = Point(); m_pocaVector = Vector();
    }

    double m_dist;        // The distance of closest approach
    double m_arcLength;   // Length along the track to the poca
    int m_region;         // One of the enums defined in Event/Recon/AcdRecon/AcdTkrPoca.h
    Point m_poca;         // Point of closest approach
    Vector m_pocaVector;  // Vector from Track to POCA
  };

 public:
  
  // retrieve Gaudi interface ID
  static const InterfaceID& interfaceID()
    { return IID_AcdIPocaTool; }
  
  AcdIPocaTool() {}
  virtual ~AcdIPocaTool() {}
    
  // @brief calculate the distance of closest approach between the track and the tile center
  virtual StatusCode doca (const AcdTileDim& tile,
			   const Event::TkrTrack& aTrack, 
			   PocaData& data) = 0;
  
  // @brief calculate the distance in the plane of the tile to the nearest edge of the tile
  //  this calculation is positive if the track hits the tile and negative otherwise
  //  the returned value is the largest number for any of the four tile edges
  virtual StatusCode hitTileDist(const AcdTileDim& tile,
				 const Event::TkrTrack& aTrack, 
				 PocaData& data) = 0;

  // @brief calculate the 3D distance to the nearest edge of the tile.  This calculation is 
  //  positive if the track hits the tile and negative otherwise
  //  the returned value is the largest number for any of the tile edges or corners
  virtual StatusCode tileActiveDist(const AcdTileDim& tile,
				    const Event::TkrTrack& aTrack, 
				    PocaData& data) = 0;
  
  // @brief calculate the 3D distance ribbon.  This calculation is 
  //  positive if the track hits the ribbon and negative otherwise
  virtual StatusCode hitRibbonDist(const AcdRibbonDim& ribbon,
				   const Event::TkrTrack& aTrack, 
				   PocaData& data) = 0;

  // @brief Make an AcdTrkPoca object, given the PocaData and the G4Propagator
  virtual StatusCode makePoca(const Event::TkrTrack& track, const PocaData& poca, const idents::AcdId& acdId,
			      IPropagator& g4PropTool, Event::AcdTkrPoca*& pocaCol) = 0;

} ;

#endif




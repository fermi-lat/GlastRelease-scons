
#ifndef __AcdTkrIntersectTool_H
#define __AcdTkrIntersectTool_H 1

#include "AcdITkrIntersectTool.h"
#include "AcdIPocaTool.h"
#include "AcdGeomMap.h"

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
* @class AcdTkrIntersectTool
*
* Class for extrapolating the Tracks to their intersections w/ the ACD
*
* $ $
*/


class AcdTkrIntersectTool : public AcdITkrIntersectTool,  public AlgTool {
	
 public:
    
  AcdTkrIntersectTool
    ( const std::string & type, 
      const std::string & name,
      const IInterface * parent ) ;
  virtual ~AcdTkrIntersectTool() ;
  
  /// @brief Intialization of the tool
  virtual StatusCode initialize();

  /// @brief main method
  virtual StatusCode makeIntersections(IPropagator& prop,
				       const AcdRecon::TrackData& track,
				       const AcdRecon::ExitData& data,	
				       const AcdRecon::PocaDataPtrMap& pocaMap,
				       const AcdRecon::AcdHitMap& hitMap,
				       AcdGeomMap& geomMap,
				       Event::AcdTkrIntersectionCol& intersections,
				       Event::AcdTkrGapPocaCol& gapPocas);

  // @brief calculate the arclength at which a track exits the tracking volume
  virtual StatusCode exitsLAT(const Event::TkrTrack& track, bool forward,
			      AcdRecon::ExitData& data);

  // @brief calculate the arclength at which a ray exits the tracking volume
  virtual StatusCode exitsLAT(const Point& x, const Vector& v, bool forward,
			      AcdRecon::ExitData& data);

  // @brief calculate the arclength at which a ray enters the tracking volume
  virtual StatusCode entersLAT(const Point& x, const Vector& v, bool forward,
			       AcdRecon::ExitData& data);

  // @brief make the TDS object that states where the track left the ACD
  virtual StatusCode makeTkrPoint(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
				  const Event::TkrTrackParams& params, Event::AcdTkrPoint*& tkrPoint );


protected:
  
  virtual StatusCode fallbackToNominal(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
				       Event::AcdTkrGapPocaCol& gapPocas);  
  
  virtual StatusCode holePoca(const AcdRecon::TrackData& track, const AcdRecon::PocaData& pocaData, const AcdTileDim& tile,
			      Event::AcdTkrGapPocaCol& gapPocas);  
  
  virtual StatusCode gapPocaRibbon(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
				   const AcdRecon::PocaData& pocaData, Event::AcdTkrGapPocaCol& gapPocas);
  
  virtual StatusCode gapPocaTile(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
				 const AcdRecon::PocaData& pocaData, AcdGeomMap& geomMap, Event::AcdTkrGapPocaCol& gapPocas);

  virtual StatusCode gapPocaCorner(const AcdRecon::TrackData& track, const AcdRecon::ExitData& data,
				   Event::AcdTkrGapPocaCol& gapPocas);

  virtual StatusCode makeGapPoca(idents::AcdGapId& gapId, const AcdRecon::TrackData& track, const AcdRecon::PocaData& pocaData,
				 double distance, Event::AcdTkrGapPoca*& poca);


  static bool checkVolId(idents::VolumeIdentifier& volId);

private:

  AcdIPocaTool*    m_pocaTool;
  IPropagator *    m_G4PropTool; 
  IGlastDetSvc*    m_detSvc; 
  IAcdGeometrySvc* m_acdGeomSvc;

  // Define the fiducial volume of the LAT
  // FIXME -- this should come for some xml reading service
  //
  // top is defined by planes at + 754.6 -> up to stacking of tiles
  // sides are defined by planes at +-840.14
  // the bottom of the ACD is at the z=-50 plane

  // Later we add 10 cm to make sure that we catch everything
  static const double s_top_distance; // = 754.6;     // center of tiles in cols 1 and 3
  static const double s_side_distance; // = 840.14;   // center of tiles in sides
  static const double s_bottom_distance; // = -50.; 
 
} ;

#endif
	
	
	

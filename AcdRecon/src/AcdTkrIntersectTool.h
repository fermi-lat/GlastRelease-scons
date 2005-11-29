
#ifndef __AcdTkrIntersectTool_H
#define __AcdTkrIntersectTool_H 1

#include "AcdITkrIntersectTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GlastSvc/Reco/IPropagatorTool.h"
#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/Reco/IPropagator.h" 
#include "TkrUtil/ITkrGeometrySvc.h"

#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrTrackParams.h"

class MsgStream;
class HepMatrix;

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

  /// @brief Default cluster finding framework
  virtual StatusCode findIntersections (const Event::TkrTrackCol *,
					Event::AcdTkrIntersectionCol *,
					std::map<idents::AcdId,unsigned char>& );
  
 protected:

  virtual int doTrack(const Event::TkrTrack& aTrack, int iTrack, std::map<idents::AcdId,unsigned char>& hitMap, MsgStream& log,
		      bool forward = true);

  void errorAtXPlane(double delta, const Event::TkrTrackParams& track, HepMatrix& covAtPlane) const;
  void errorAtYPlane(double delta, const Event::TkrTrackParams& track, HepMatrix& covAtPlane) const;
  void errorAtZPlane(double delta, const Event::TkrTrackParams& track, HepMatrix& covAtPlane) const;

 private:

  Event::AcdTkrIntersectionCol* output;

  IPropagator *    m_G4PropTool; 
  IGlastDetSvc*    m_detSvc; 
  ITkrGeometrySvc* m_tkrGeom;

} ;

#endif
	
	
	

#ifndef __ACD_RECONV2_H
#define __ACD_RECONV2_H 1
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/Digi/AcdDigi.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"
#include "Event/Recon/AcdRecon/AcdReconV2.h"
#include "Event/Recon/AcdRecon/AcdEventTopology.h"
#include "Event/Recon/AcdRecon/AcdPocaMap.h"
#include "Event/Recon/AcdRecon/AcdTkrHitPoca.h"
#include "Event/Recon/AcdRecon/AcdTkrGapPoca.h"
#include "Event/Recon/AcdRecon/AcdHit.h"
#include "Event/Recon/AcdRecon/AcdAssoc.h"

#include "TkrUtil/ITkrTrackVecTool.h"

#include "GaudiKernel/ObjectVector.h"

#include "AcdITkrIntersectToolV2.h"
#include "AcdIPha2MipTool.h"
#include "AcdIPocaToolV2.h"

#include "AcdUtil/AcdGeomMap.h"

#include "../AcdRecon/AcdReconStruct.h"
#include "../AcdRecon/AcdPatRecTools.h"

#include "GlastSvc/Reco/IPropagatorTool.h"
#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/Reco/IPropagator.h" 
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "AcdUtil/IAcdGeometrySvc.h"
#include "idents/AcdId.h"
#include "idents/VolumeIdentifier.h"
#include "geometry/Ray.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
// TU: Hacks for CLHEP 1.9.2.2 and beyond
#ifndef HepPoint3D
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#ifndef HepVector3D
typedef HepGeom::Vector3D<double> HepVector3D;
#endif

#include <map>
#include <vector>

/** @class AcdReconAlgV2
 * @brief ACD reconstruction using the AcdDigi collection from the TDS.
 *
 * 
 *
 * TDS Inputs:
 *  - Event::AcdDigiCol: all the AcdDigi objects
 *  - Event::McParticleCol: MC tracks, obviously for MC only
 *
 * TDS Outputs:
 *  - Event::AcdRecon object, which contains
 *    - Event::AcdTkrIntersectionCol: all the intersections with GEANT model 
 *    - Event::AcdTkrHitPocaCol all the POCA calculations
 *    - Event::AcdTkrGapPocaCol: POCA w.r.t. gaps and ribbons in the GEANT model 
 *    - Event::AcdTkrPointCol: Data about where the tracks exits the ACD volume
 *    - Event::AcdSplashVarsCol: NOT Filled!! (Data about backsplash from tracks)
 *    - Some numbers we extracted during the recon process
 *      - Number of hit tiles
 *      - Number of hit ribbons
 *      - Max active distance over all tiles, ID of tile w/ max active distance
 *      - Max active distance over all ribbons, ID of ribbon w/ max active distance
 *      - Max active distance for all rows of tiles of ACD
 *      - Minimum distance to a corner ray
 *    - Some MC quantaties
 *      - Total tile energy (MC)
 *      - Total ribbon energy (MC)
 *      - twin vectors of AcdID and MC energy
 *    - Some useless crap
 *      - Gamma DOCA (Deprecated!!)
 *      - Min DOCA over all tracks, ID of tile w/ min DOCA (Deprecated!!)
 *      - Min DOCA and ID for all rows of ACD (Deprecated!!)
 *
 *
 * Note: In some collections the objects may differ in substantial ways.
 *
 *  - Track extensions can go in either direction:
 *     - arcLength > 0 -> up
 *     - arcLength < 0 -> down
 *  - Objects have index of underlying track:
 *     - track_index >= 0 -> track
 *     - track_index = -1 -> vertex
 * 
 * 
 * Algorithm::
 *  - For all data:
 *    -# Use AcdPha2MipTool to make fill the Event::AcdHitCol object
 *    -# Fill MC energy and tile ID vectors
 *    -# Calculate all the track related distances in trackDistances()
 *      - Loop on tracks
 *        - Extend tracks in both directions, find LAT exit point using AcdRecon::exitsLat()
 *        - Use AcdPocaTool to calculate POCA w.r.t. all hits in Event::AcdHitCol in both directions
 *        - Loop on POCA objects, latch largest active distances in tileActDist() and hitRibbonDist()
 *          - do tiles in both directions, ribbons downwards only
 *        - Use AcdPocaTool to filter POCA by activeDistance
 *        - Use AcdTkrIntersectionTool to extrapolate track in both directions, fills:
 *          - Event::AcdTkrIntersectionCol
 *          - Event::AcdPocaSet
 *          - Event::AcdTkrGapPocaCol
 *          - Event::AcdTkrPointCol
 *          
 *    -# Calculate all the vertex related distances
 *      - Extend vertex fit in both directions, find LAT exit point using AcdRecon::exitsLat()
 *      - Use AcdPocaTool to calculate POCA w.r.t. all hits in Event::AcdHitCol in both directions
 *      - Use AcdPocaTool to filter POCA by activeDistance
 *      - Use AcdTkrIntersectionTool to extrapolate vertex in both directions, fills:
 *        - Event::AcdPocaSet
 *        - Event::AcdTkrPointCol
 *
 *    -# Sort all the AcdRecon::PocaData objects ( largest active distance first )
 *    -# Use AcdPocaTool to make Event::AcdTkrPocaData objects
 *    -# Calculate backsplash variable (NO-OP)
 *    -# Put data on TDS under Event/AcdRecon
 *
 *  - For MC only:
 *    -# Calculate all the MC track related distances
 *      - Extend MC track through LAT
 *      - Use AcdPocaTool to calculate POCA w.r.t. all hits in Event::AcdHitCol in both directions
 *      - Use AcdPocaTool to filter POCA by activeDistance
 *      - Use AcdTkrIntersectionTool to extrapolate vertex upward to MC entry point
 *    -# Sort all the MC related AcdRecon::PocaData objects ( largest active distance first )
 *    -# Use AcdPocaTool to make Event::AcdTkrPocaData objects
 *    -# Put data on TDS under Event/MC/McAcdTkrHitPocaCol Event/MC/McAcdTkrPointCol
 *
 * 
 * This Algorithm has 5 JO paramters:  
 *  - Tool names:
 *    -  intersectionToolName["AcdTkrIntersectTool"]
 *    -  hitToolName["AcdPha2MipTool"]
 *    -  pocaToolName["AcdPocaTool"]
 *    -  propToolName["G4PropagationTool"]
 *  - doBackSplash[false] : do turn on backsplash calculations (caveat emptor)
 *
 * @author Heather Kelly
 * @author Eric Charles
 *
 * $Header$
 */
class AcdReconAlgV2 : public Algorithm
{

  public:
      AcdReconAlgV2(const std::string& name, ISvcLocator* pSvcLocator); 

      StatusCode initialize();
      StatusCode execute();
      StatusCode finalize();
  
  private:

      struct hit_pointer_less : public std::binary_function<Event::AcdTkrHitPoca*,Event::AcdTkrHitPoca*, bool> {
	bool operator()(Event::AcdTkrHitPoca* x, Event::AcdTkrHitPoca* y) { 	  
	  return x->operator<(*y);
	}
      };
      
      struct gap_pointer_less : public std::binary_function<Event::AcdTkrGapPoca*,Event::AcdTkrGapPoca*, bool> {
	bool operator()(Event::AcdTkrGapPoca* x, Event::AcdTkrGapPoca* y) { 
	  return x->operator<(*y); }
      };

      /// reset all member variables for each iteration
      void clear ();

      /// Retrieve geometry parameters
      void getParameters ();
      
      /// fill the pattern recognition map
      StatusCode fillPatRecMap( );

      /// routine called by execute that performs the reconstruction 
      StatusCode reconstruct (const Event::AcdDigiCol& digiCol);

      /// routine called by execute that performs reconstruction on MC side
      StatusCode doMC(const Event::AcdHitCol& acdHits);

      /// retrieves MC particles and calls the DOCA and Active Distance routines
      StatusCode mcDistances(const Event::AcdHitCol& acdHits, 
			     Event::AcdTkrAssocCol& tkrAssocs);

      /// retrieves tracks and calls the DOCA and Active Distance routines
      StatusCode trackDistances(const Event::AcdHitCol& acdHits, 
				Event::AcdTkrAssocCol& tkrAssocs,
				bool upward = true);

      /// try and do the same thing with the CAL clusters
      StatusCode calClusterDistances(const Event::AcdHitCol& acdHits, 
				     Event::AcdCalAssocCol& calAssocs);

      /// retrieves event vertex and calls the DOCA and Active Distance routines
      StatusCode vertexDistances(const Event::AcdHitCol& acdHits, 
				 Event::AcdTkrAssocCol& tkrAssocs,
				 bool upward = true);

      /// Calculates the point where the paritcle crosses the nominal ACD 
      StatusCode exitPoint(const AcdRecon::TrackData& aTrack, bool forward,
			   AcdRecon::ExitData& data, double tolerance = 0.);

      /// get the all the distances to hit tiles & ribbons for track in one direction
      StatusCode hitDistances(const AcdRecon::TrackData& aTrack, const Event::AcdHitCol& acdHits, 
			      const AcdRecon::ExitData& data,
			      AcdRecon::PocaDataMap& pocaMap);				  

      /// get the distances for a single element and a single track in one detector
      StatusCode elemDistances(const AcdRecon::TrackData& aTrack, const idents::AcdId& acdId,
			       AcdRecon::PocaData& pocaData);

      /// Extrapolate track as far as needed, add error to AcdTkrPoca, make AcdTkrIntersections
      StatusCode extrapolateTrack(const Event::TkrTrackParams& trackParams,
				  const AcdRecon::TrackData& trackData,
				  const AcdRecon::ExitData& isectData,
				  AcdRecon::PocaDataPtrMap& pocaDataMap,
				  int& ssdVeto,
				  std::vector<Event::AcdTkrHitPoca*>& hitPocae,
				  std::vector<Event::AcdTkrGapPoca*>& gapPocae,
				  Event::AcdTkrPoint*& point);

      ///  Extrapolate & build TDS object
      StatusCode extrapolateVertex(const AcdRecon::TrackData& trackData,
				   const AcdRecon::ExitData& isectData,
				   AcdRecon::PocaDataPtrMap& pocaDataMap,
				   int& ssdVeto,
				   std::vector<Event::AcdTkrHitPoca*>& hitPocae,
				   std::vector<Event::AcdTkrGapPoca*>& gapPocae,
				   Event::AcdTkrPoint*& points);

      /// Pick the best pocae for a track
      StatusCode sortPocae(std::vector<Event::AcdTkrHitPoca*>& hitPocae, 
			   std::vector<Event::AcdTkrGapPoca*>& gapPocae );

      /// Fill an AcdTkrAssoc with data
      StatusCode fillTkrAssoc(Event::AcdAssoc& assoc,
			      const std::vector<Event::AcdTkrHitPoca*>& hitPocae,
			      const std::vector<Event::AcdTkrGapPoca*>& gapPocae,
			      Event::AcdTkrPoint* point);

      /// Fill an AcdCalAssoc with data
      StatusCode fillCalAssoc(Event::AcdAssoc& assoc,
			      const std::vector<Event::AcdTkrHitPoca*>& hitPocae,
			      const std::vector<Event::AcdTkrGapPoca*>& gapPocae,
			      Event::AcdTkrPoint* point);

      /// Fill the AcdEventTopology object
      StatusCode fillAcdEventTopology(const Event::AcdHitCol& acdHits,
				      Event::AcdEventTopology&);
				      

      StatusCode calcCornerDoca(const AcdRecon::TrackData& trackData,
                                float &dist);

      /// the tool to calculate the Track intersections w/ the ACD
      AcdITkrIntersectToolV2* m_intersectionTool;

      /// name of Tool for finding the Track intersections
      std::string m_intersectionToolName;

      /// the tool to calculate the AcdHits in terms on MIPS
      AcdIPha2MipTool* m_hitTool;

      /// name of Tool for makeint the hits
      std::string m_hitToolName;

      /// the tool to calculate the AcdTkrPocas
      AcdIPocaToolV2* m_pocaTool;
      
      /// name of Tool for makeint the hits
      std::string m_pocaToolName;
      
      IPropagator *    m_G4PropTool; 

      std::string m_propToolName;

      /// variables to store instrument parameters
      static double s_vetoThresholdMeV;
      static unsigned int s_numSideRows;
      static AcdRecon::AcdVolume s_acdVolume;

      /// access to the Glast Detector Service to read in geometry constants from XML files
      IGlastDetSvc *m_glastDetSvc;
      IAcdGeometrySvc *m_acdGeoSvc;

      /// trackVecTool
      ITkrTrackVecTool* m_pTrackVec;
  
      /// Tolernace for pat-rec hash map
      double           m_patRecTol;

      /// Turn on or off the calculation of the Backsplash variables
      bool m_doBackSplash;

      /// Turn on or off the extrapolations using the vertex solutions
      bool m_doVertexExtrap;
      
      /// Turn on or off the extrapolations using the downgoing solutions
      bool m_doDownwardExtrap;
      
      /// map of AcdId and the corresponding hit status 
      AcdRecon::AcdHitMap m_hitMap;

      /// map to keep track to the tile and ribbon geoms
      AcdGeomMap* m_geomMap;

      /// The map used for pattern recognition
      AcdRecon::AcdPatRecMap* m_patRecMap;



};

#endif

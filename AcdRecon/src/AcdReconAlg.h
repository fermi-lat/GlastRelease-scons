#ifndef __ACD_RECON_H
#define __ACD_RECON_H 1
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/Digi/AcdDigi.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"
#include "Event/Recon/AcdRecon/AcdPocaMap.h"
#include "Event/Recon/AcdRecon/AcdTkrGapPoca.h"
#include "Event/Recon/AcdRecon/AcdHit.h"

#include "GaudiKernel/ObjectVector.h"

#include "AcdITkrIntersectTool.h"
#include "AcdIPha2MipTool.h"
#include "AcdIPocaTool.h"

#include "AcdUtil/AcdGeomMap.h"

#include "../AcdRecon/AcdReconStruct.h"

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

/** @class AcdReconAlg
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
class AcdReconAlg : public Algorithm
{

  public:
      AcdReconAlg(const std::string& name, ISvcLocator* pSvcLocator); 

      StatusCode initialize();
      StatusCode execute();
      StatusCode finalize();
  
  private:

      /// reset all member variables for each iteration
      void clear ();

      /// Retrieve geometry parameters
      void getParameters ();

      /// routine called by execute that performs the reconstruction 
      StatusCode reconstruct (const Event::AcdDigiCol& digiCol);

      /// routine called by execute that performs reconstruction on MC side
      StatusCode doMC(const Event::AcdHitCol& acdHits);

      /// retrieves MC particles and calls the DOCA and Active Distance routines
      StatusCode mcDistances(const Event::AcdHitCol& acdHits, 
			     Event::AcdPocaSet& pocaSet,
			     Event::AcdTkrPointCol& exitPoints);

      /// retrieves tracks and calls the DOCA and Active Distance routines
      StatusCode trackDistances(const Event::AcdHitCol& acdHits, 
				Event::AcdPocaSet& pocaSet,
				Event::AcdTkrIntersectionCol& acdIntersections,
				Event::AcdTkrGapPocaCol& gapPocas,
				Event::AcdTkrPointCol& exitPoints);

      /// retrieves event vertex and calls the DOCA and Active Distance routines
      StatusCode vertexDistances(const Event::AcdHitCol& acdHits, 
				 Event::AcdPocaSet& pocaSet,
				 Event::AcdTkrPointCol& exitPoints);

      /// Calculates the point where the paritcle crosses the nominal ACD 
      StatusCode exitPoint(const AcdRecon::TrackData& aTrack, bool forward,
			   AcdRecon::ExitData& data, double tolerance = 0.);

      /// get the all the distances to hit tiles & ribbons for track in one direction
      StatusCode hitDistances(const AcdRecon::TrackData& aTrack, const Event::AcdHitCol& acdHits, 
			      AcdRecon::PocaDataMap& pocaMap);				  

      /// Bill Atwood's new calculation for Active Distance, in 3D
      StatusCode tileActDist(const AcdRecon::PocaDataMap& pocaMap,
			     std::vector<double> &row_values, double &dist, 
                             idents::AcdId& maxActDistId);

      /// Bill Atwood's new calculation for Active Distance - applied to ribbons
      StatusCode hitRibbonDist(const AcdRecon::PocaDataMap& pocaMap,
			       double &dist, idents::AcdId& maxActDistId);

      /// Extrapolate track as far as needed, add error to AcdTkrPoca, make AcdTkrIntersections
      StatusCode extrapolateTrack(const Event::TkrTrack& aTrack,
				  const AcdRecon::TrackData& trackData,
				  const AcdRecon::PocaDataPtrMap& pocaDataMap,
				  const AcdRecon::ExitData& isectData,
				  Event::AcdPocaSet& pocaSet,
				  Event::AcdTkrIntersectionCol& acdIntersections,
				  Event::AcdTkrGapPocaCol& gapPocas,
				  Event::AcdTkrPointCol& points);

      ///  Extrapolate & build TDS object
      StatusCode extrapolateVertex(const AcdRecon::TrackData& trackData,
				   const AcdRecon::PocaDataPtrMap& pocaDataMap,
				   const AcdRecon::ExitData& isectData,
				   Event::AcdPocaSet& pocaSet,
				   Event::AcdTkrPointCol& points);

      /// Write an exit point object
      void writeExitPoint(std::ostream& os, 
			  const AcdRecon::TrackData& trackData, const AcdRecon::ExitData& isectData);			  

      StatusCode calcCornerDoca(const HepPoint3D &x0, const HepVector3D &dir,
                                double &dist);

      /// calculate stuff for the backsplash
      StatusCode doBacksplash(const Event::AcdDigiCol& digiCol, Event::AcdSplashVarsCol& acdSplashVars);

      /// the tool to calculate the Track intersections w/ the ACD
      AcdITkrIntersectTool* m_intersectionTool;

      /// name of Tool for finding the Track intersections
      std::string m_intersectionToolName;

      /// the tool to calculate the AcdHits in terms on MIPS
      AcdIPha2MipTool* m_hitTool;

      /// name of Tool for makeint the hits
      std::string m_hitToolName;

      /// the tool to calculate the AcdTkrPocas
      AcdIPocaTool* m_pocaTool;
      
      /// name of Tool for makeint the hits
      std::string m_pocaToolName;
      
      IPropagator *    m_G4PropTool; 

      std::string m_propToolName;

      /// variables to store instrument parameters
      static double s_vetoThresholdMeV;
      static unsigned int s_numSideRows;
      static AcdRecon::AcdVolume s_acdVolume;

      // record of the tile with the minimum Distance of Closest Approach
      idents::AcdId m_minDocaId, m_ribbon_act_dist_id, m_maxActDistId,  
                    m_maxActDist3DId, m_maxActDist3DId_down;

      /// access to the Glast Detector Service to read in geometry constants from XML files
      IGlastDetSvc *m_glastDetSvc;
      IAcdGeometrySvc *m_acdGeoSvc;
  
    
      bool m_calcCornerDoca;

      bool m_doBackSplash;

      /// Number of Acd Tiles above threshold
      unsigned int m_tileCount, m_ribbonCount;
      /// Total Energy deposited in the ACD system
      double m_totEnergy, m_gammaDoca, m_totRibbonEnergy;
      /// Minimun Distance of Closest Approach
      double m_doca;
      /// DOCA to corner gaps
      double m_cornerDoca;
      /// Minimum Active Distance
      double m_act_dist, m_ribbon_act_dist, m_act_dist3D, m_act_dist3D_down;
      /// list of DOCA values for top and each side row
      std::vector<double> m_rowDocaCol;
      /// list of active distance values for top and each side row
      std::vector<double> m_rowActDistCol;
      std::vector<double> m_rowActDist3DCol;
      std::vector<double> m_rowActDist3DCol_down;
      /// map of AcdId and their corresponding energies
      //std::map<idents::AcdId, double> m_energyCol;
      std::vector<idents::AcdId> m_idCol, m_idRibbonCol;
      std::vector<double> m_energyCol, m_energyRibbonCol;
   
      /// map of AcdId and the corresponding hit status 
      AcdRecon::AcdHitMap m_hitMap;

      /// map to keep track to the tile and ribbon geoms
      AcdGeomMap* m_geomMap;

};

#endif

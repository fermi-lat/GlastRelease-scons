#ifndef __ACD_RECON_H
#define __ACD_RECON_H 1
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/Digi/AcdDigi.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"
#include "Event/Recon/AcdRecon/AcdPocaMap.h"
#include "Event/Recon/AcdRecon/AcdTkrGapPoca.h"

#include "GaudiKernel/ObjectVector.h"

#include "AcdITkrIntersectTool.h"
#include "AcdIPha2MipTool.h"
#include "AcdIPocaTool.h"

#include "AcdGeomMap.h"

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
 * Computes a number of quantities that are then stored in the AcdRecon object
 * and put on the TDS.  Those quantities that are computed includes:
 * - Minimum Distance of Closest Approach (DOCA)
 * - List of DOCA values containing the min DOCA for each row and top ACD tiles.
 * - Minimum Active Distance quantity for all hit ACD tiles
 * - List of Active Distance values containing min. Active Distance for each 
 * row and the top tiles.
 *
 * The DOCA and Active Distance quantities are computed using the ACD detector hits
 * and the TkrRecon reconstructed track collection.  DOCA is calculated by finding
 * the minimum distance between the center of hit ACD tiles and all found tracks.
 * Active Distance is calculated by finding the minimum distance between the edge
 * of hit ACD tiles and all found tracks.
 *
 * @author Heather Kelly
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

      /// retrieves tracks and calls the DOCA and Active Distance routines
      StatusCode trackDistances(const Event::AcdDigiCol& digiCol, 
				Event::AcdPocaSet& pocaSet,
				Event::AcdTkrIntersectionCol& acdIntersections,
				Event::AcdTkrGapPocaCol& gapPocas,
				Event::AcdTkrPointCol& exitPoints);

      /// Calculates the point where the paritcle crosses the nominal ACD 
      StatusCode exitPoint(const AcdRecon::TrackData& aTrack, bool forward,
			   AcdRecon::ExitData& data, double tolerance = 0.);

      /// get the all the distances to hit tiles & ribbons for track in one direction
      StatusCode hitDistances(const AcdRecon::TrackData& aTrack, const Event::AcdDigiCol& digiCol, 
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

      // record of the tile with the minimum Distance of Closest Approach
      idents::AcdId m_minDocaId, m_ribbon_act_dist_id, m_maxActDistId,  
                    m_maxActDist3DId;

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
      double m_act_dist, m_ribbon_act_dist, m_act_dist3D;
      /// list of DOCA values for top and each side row
      std::vector<double> m_rowDocaCol;
      /// list of active distance values for top and each side row
      std::vector<double> m_rowActDistCol;
      std::vector<double> m_rowActDist3DCol;
      /// map of AcdId and their corresponding energies
      //std::map<idents::AcdId, double> m_energyCol;
      std::vector<idents::AcdId> m_idCol, m_idRibbonCol;
      std::vector<double> m_energyCol, m_energyRibbonCol;
   
      /// map of AcdId and the corresponding hit status 
      AcdRecon::AcdHitMap m_hitMap;

      /// map to keep track to the tile and ribbon geoms
      AcdGeomMap m_geomMap;

};

#endif

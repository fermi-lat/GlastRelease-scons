#ifndef __ACD_RECON_H
#define __ACD_RECON_H 1
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/Digi/AcdDigi.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"

#include "GaudiKernel/ObjectVector.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "idents/AcdId.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"

#include <map>

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
      StatusCode trackDistances();

      /// Old style - distance of closest approach calculation
      /// Finds minimum perpendicular distance from tracks to the center of the tiles
      double doca (const HepPoint3D &x0, const HepVector3D &dir, std::vector<double> &doca_values);

      /// Bill Atwood's new calculation for Active Distance
      double hitTileDist(const HepPoint3D &x0, const HepVector3D &dir, std::vector<double> &row_values);

      /// variables to store instrument parameters
      static double s_vetoThresholdMeV;
      static unsigned int s_numSideRows;

      // record of the tile with the minimum Distance of Closest Approach
      idents::AcdId m_minDocaId;

      /// access to the Glast Detector Service to read in geometry constants from XML files
      IGlastDetSvc *m_glastDetSvc;

      Event::AcdDigiCol m_acdDigiCol;

      Event::AcdRecon *m_acdRecon;

      /// Number of Acd Tiles above threshold
      unsigned int m_tileCount;
      /// Total Energy deposited in the ACD system
      double m_totEnergy, m_gammaDoca;
      /// Minimun Distance of Closest Approach
      double m_doca;
      /// Minimum Active Distance
      double m_act_dist;
      /// list of DOCA values for top and each side row
      std::vector<double> m_rowDocaCol;
      /// list of active distance values for top and each side row
      std::vector<double> m_rowActDistCol;
      /// map of AcdId and their corresponding energies
      //std::map<idents::AcdId, double> m_energyCol;
	  std::vector<idents::AcdId> m_idCol;
	  std::vector<double> m_energyCol;
};

#endif

#ifndef __ACD_RECON_H
#define __ACD_RECON_H 1
#include "GaudiKernel/Algorithm.h"

// Glast specific includes
#include "Event/TopLevel/Event.h"
#include "Event/Digi/AcdDigi.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/AcdRecon.h"

#include "GaudiKernel/ObjectVector.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "idents/AcdId.h"
#include "geometry/Vector.h"
#include <map>

class Point;
class Ray;

/** @class AcdReconAlg
 * @brief ACD reconstruction
 *
 * Migration of old VetoRecon code in glastsim to a Gaudi algorithm.
 * DOCA stands for Distance of Closest Approach
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

      /// retrieves tracks and calls the DOCA routines
      StatusCode acdDoca();

      /// Old style - distance of closest approach calculation
      /// Finds minimum perpendicular distance from tracks to the center of the tiles
      double doca (const Point &x0, const Vector &dir, std::vector<double> &doca_values);

      /// Bill Atwood's new calculation for Active Distance
      double hitTileDist(const Point &x0, Vector &dir);

      /// variables to store instrument parameters
      static double s_thresholdEnergy;
      static unsigned int s_numSideRows;

      // record of the tile with the minimum Distance of Closest Approach
      idents::AcdId m_minDocaId;

      /// access to the Glast Detector Service to read in geometry constants from XML files
      IGlastDetSvc *m_glastDetSvc;

      Event::AcdDigiCol m_acdDigiCol;

      Event::AcdRecon *m_acdRecon;

      /// Items that will be output to the ntuple
      unsigned int m_tileCount;
      double m_totEnergy, m_gammaDoca, m_doca, m_act_dist;
      std::vector<double> m_rowDocaCol;
      std::map<idents::AcdId, double> m_energyCol;


};

#endif

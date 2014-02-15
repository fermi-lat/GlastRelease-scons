#ifndef _AcdDigi_AcdDigiAlg_H
#define _AcdDigi_AcdDigiAlg_H 1

#include "GaudiKernel/Algorithm.h"

#include "idents/AcdId.h"

#include "AcdTileList.h"
#include "AcdDigiUtil.h"

#include <map>

#include "Event/Digi/AcdDigi.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "GlastSvc/GlastRandomSvc/IRandomAccess.h"


#include "AcdUtil/IAcdFailureModeSvc.h"

/** @class AcdDigiAlg
* @brief Algorithm to convert from hit data stored as McPositionHits into digitization data 
* for the ACD.
* 
* @author Heather Kelly
*
* $Header$
*/

class AcdDigiAlg : public Algorithm {
    
public:
    
  AcdDigiAlg(const std::string& name, ISvcLocator* pSvcLocator); 
  
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

protected:

  StatusCode fillEnergyAndPeMaps( const Event::McPositionHitCol& mcHits,
                 std::map<idents::AcdId, double>& energyIdMap,
                 std::map<idents::AcdId, std::pair<double, double> >& peMap,
                 std::map<idents::AcdId, unsigned int>& statusMap);
  
  StatusCode convertPeToMips( const std::map<idents::AcdId, std::pair<double, double> >& peMap,
                              const std::map<idents::AcdId, unsigned int>&               statusMap,
                                    std::map<idents::AcdId, std::pair<double, double> >& mipsMap);
  
  StatusCode makeDigis(const std::map<idents::AcdId, std::pair<double,double> >& mipsMap,
                       const std::map<idents::AcdId, unsigned int>& statusMap,
                       Event::AcdDigiCol& digiCol);
  
  /// Clear the local paramaters
  void clear();  
        
private:

  IAcdFailureModeSvc *m_acdFailureSvc;

  IRandomAccess *m_randTool;

  /// input XML file containing parameters for Digitization
  std::string m_xmlFileName;
        
  /// input AcdCalibration service
  std::string m_calibSvcName;  

  /// JobOptions parameter denoting whether or not to apply Poisson fluctuations
  bool m_apply_poisson;
  
  /// JobOptions parameter denoting whether or not to apply Gaussian noise 
  /// before determining PHA and discriminators
  bool m_apply_noise;
  
  /// JobOptions parameter denoting whether or not to apply Coherent readout noise 
  /// to the PHA values
  bool m_apply_coherent_noise;

  /// JobOptions parameter denoting whether or not coherent noise is in the overlay events 
  /// to the PHA values
  bool m_coherent_noise_in_overlay;

  /// JobOptions parameter denoting whether or not to apply edge effects
  /// according to the position of MC hits.
  bool m_edge_effect;

  /// JobOptions parameter use to degrade performance of ACD
  /// Reduces the light yield, defaults to 1.0
  double m_lightYeildRatio;

  /// JobOptions paramter denoting zero level suppression
  double m_phaZeroThreshold;

  /// Access the methods in the AcdDigiUtil class
  AcdDigiUtil m_util;
  
  /// A list of all the tiles and ribbons
  AcdTileList m_tiles;
    
  std::map<idents::AcdId, double>                     m_energyDepMap;
  std::map<idents::AcdId, std::pair<double, double> > m_peMap;
  std::map<idents::AcdId, std::pair<double, double> > m_mipsMap;
  std::map<idents::AcdId, unsigned int>               m_statusMap;

  unsigned short m_gemDeltaEventTime;


};

#endif

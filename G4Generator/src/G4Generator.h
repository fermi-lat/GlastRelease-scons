#ifndef G4Generator_h
#define G4Generator_h

#include "GaudiKernel/Algorithm.h"

#include <string>
#include <vector>
#include "GaudiKernel/Property.h"

class RunManager;
class IParticlePropertySvc;
class IG4GenErrorSvc;


/** 
 * @class G4Generator
 *
 * @brief Gaudi algorithm for GLAST simulation
 *
 * This is the main Gaudi algorithm to use the Geant4 toolkit as a Montecarlo
 * propagator
 *  
 * @author T.Burnett and R.Giannitrapani
 *    
 * $Header$
 */
class G4Generator : public Algorithm {
 public:
  G4Generator(const std::string& name, ISvcLocator* pSvcLocator); 
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();
 
  
 private:
  /// a pointer to the service for particle properties
  IParticlePropertySvc* m_ppsvc;

  /// a pointer to the Error handling service
  IG4GenErrorSvc* m_ErrorSvc;

  /// set of UI commands for setup; this is a property of the algorithm and can
  /// be set in the jobOptions file
  StringArrayProperty m_uiCommands;
  
  /// This is the G4 manager that handles the simulation
  RunManager* m_runManager;

  /// the geometry level of details
  std::string m_geometryMode;
  
  /// flag to save trajectories in the TDS
  bool m_saveTrajectories;

  /// Maximum number of generations to store trajectories for
  int m_numGenerations;
  double m_minDistance;
  double m_lowEnergy;

  /// the McParticle tree mode
  /// It can be "full" or "minimal" or "pruneCal"

  std::string m_mcTreeMode;
  
  /// The default cutoff value (in mm)
  /// This now includes cutoffs for the tracker, calorimeter and 
  /// everything else
  DoubleProperty m_defaultCutValue;
  DoubleProperty m_defaultTkrCutValue;
  DoubleProperty m_defaultCalCutValue;

  /// the Physics List
  /// It can be "full" or "only_em"

  std::string m_physics_choice;

  /// the Physics Tables could be built at initialisation time ("build")
  /// or stored in the appropriate directory ("store") or retrieved from
  /// the specified directory ("retrieve") 
  /// The directory name shall be specified by the user

  std::string m_physics_table;
  std::string m_physics_dir;

  /// if true, use Glast special version
  BooleanProperty m_mscatOption;
  BooleanProperty m_eLossCurrent;

  /// for material printout
  BooleanProperty m_printRadLen;

  /// cut in energy
  
  DoubleProperty m_defaultHeCut;

};

#endif







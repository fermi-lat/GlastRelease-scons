#ifndef G4Generator_h
#define G4Generator_h

#include "GaudiKernel/Algorithm.h"

#include <string>
#include <vector>
class IFluxSvc;
class IFlux;
class IParticlePropertySvc;
class RunManager;
namespace gui{class GuiMgr;}

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
 
  /// This method is used to setup the link with the FluxSvc to retrive the
  /// primary particle 
  /// @param source_name The name of the source in the XML sources 
  /// definition file
  void setSource(std::string source_name);
  
 private:
  /// a pointer to the flux service main classes
  IFluxSvc* m_fluxSvc;
  IFlux*    m_flux;

  IParticlePropertySvc* m_ppsvc;

  /// source name to get from the Flux service
  std::string m_source_name;

  /// set of UI commands for setup
  StringArrayProperty m_uiCommands;
  
  /// This is the G4 manager that handles the simulation
  RunManager* m_runManager;

  /// access to the GuiManager
  gui::GuiMgr* m_guiMgr;
  void setupGui();
};

#endif





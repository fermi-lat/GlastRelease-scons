// $Header$


#ifndef G4Generator_h
#define G4Generator_h

#include "GaudiKernel/Algorithm.h"

#include <string>
#include <vector>
class IFluxSvc;
class IFlux;
class RunManager;
namespace gui{class GuiMgr;}
/**
  Geant4 interface for GLAST simulation
  */

class G4Generator : public Algorithm {
public:
    G4Generator(const std::string& name, ISvcLocator* pSvcLocator); 
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
    void setSource(std::string source_name);

private:
    IFluxSvc* m_fluxSvc;
    IFlux*    m_flux;

    /// source name to get from the Flux service
    std::string m_source_name;

    /// set of UI commands for setup
    StringArrayProperty m_UIcommands;

    /// This is the G4 manager that handles the simulation
    RunManager* m_runManager;

    /// access to the GuiManager
    gui::GuiMgr* m_guiMgr;
    void setupGui();
};

#endif




// $Header$


#ifndef G4Generator_h
#define G4Generator_h

#include "GaudiKernel/Algorithm.h"

#include <string>
class IFluxSvc;
class IFlux;
class RunManager;

/**
  Geant4 interface
  */

class G4Generator : public Algorithm {
public:
    G4Generator(const std::string& name, ISvcLocator* pSvcLocator); 
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    IFluxSvc* m_fluxSvc;
    IFlux*    m_flux;

    std::string m_source_name;

    /// This is the G4 manager that handles the simulation
    RunManager* m_runManager;
};

#endif




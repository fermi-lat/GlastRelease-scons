// $Header$
#ifndef FluxAlg_h
#define FluxAlg_h

// Include files
// Gaudi system includes
#include "GaudiKernel/Algorithm.h"

class IFlux;
class IFluxSvc;
class IparticlePropertySvc;
class McVertex;

//------------------------------------------------------------------------------
/** 



*/
class FluxAlg : public Algorithm {
public:
    FluxAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
    
private: 
    
    StringProperty m_source_name;
    IFluxSvc*   m_fluxSvc;
    IFlux *     m_flux;
    
    
    unsigned long m_run;      // run number
    unsigned long m_event;    // event number
    
    McVertex*       m_root; // 
    
    IDataProviderSvc* m_eds;
    mc::McParticleCol* m_plist;

    IParticlePropertySvc * m_partSvc;
    
};
//------------------------------------------------------------------------
#endif

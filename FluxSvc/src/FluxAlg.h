// $Header$
#ifndef FluxAlg_h
#define FluxAlg_h
/** 
* \class FluxAlg
*
* \brief This is an Algorithm designed to get particle information 
* from FluxSvc and put it onto the TDS for later retrieval
* \author Toby Burnett
* 
* $Header $
*/



// Include files
// Gaudi system includes
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Property.h"

class IFlux;
class IFluxSvc;
class IparticlePropertySvc;


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
    
    
    UnsignedIntegerProperty m_run;      // run number
    unsigned int m_sequence;  // sequence number
    
    
    IDataProviderSvc* m_eds;
    
    IParticlePropertySvc * m_partSvc;
    
};
//------------------------------------------------------------------------
#endif

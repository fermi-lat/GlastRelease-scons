// $Header$
#ifndef ExposureAlg_h
#define ExposureAlg_h
/** 
* \class ExposureAlg
*
* \brief This is an Algorithm designed to get information about LAT
* position, exposure and livetime from FluxSvc and use it to put information onto the TDS about
* LAT pointing and location characteristics, effectively generating the D2 database.  The "TimeCandle" 
* Spectrum is included (and can be used in jobOptions with this algorithm) in order to provide a constant time reference.
*
* \author Sean Robinson
* 
* $Header $
*/


// Include files
// Gaudi system includes
#include "GaudiKernel/Algorithm.h"

//for file handling
#include <iostream>
#include <string>

class IFlux;
class IFluxSvc;
class IparticlePropertySvc;
class McVertex;


class ExposureAlg : public Algorithm {
public:
    ExposureAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    //stuff that an Algorithm needs.
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private: 
    
    double m_lasttime; //time value to hold time between events;
    StringProperty m_source_name;
    IFluxSvc*   m_fluxSvc;
    IFlux *     m_flux;

    //! Create the TimeCandle Spectrum.
    //void ExposureAlg::makeTimeCandle(IFluxSvc* fsvc);
    
    unsigned long m_run;      // run number
    unsigned long m_event;    // event number
    std::ostream* m_out;  //for output that looks like the stuff from the astro orbit model test.
    
    
    IDataProviderSvc* m_eds;
    
    IParticlePropertySvc * m_partSvc;
    
};
//------------------------------------------------------------------------
#endif

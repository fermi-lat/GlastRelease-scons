// $Header$
// 
//  Original author: Toby Burnett tburnett@u.washington.edu

#ifndef _H_FluxSvc_
#define _H_FluxSvc_

// includes
#include "GaudiKernel/Service.h"
#include "FluxSvc/IFluxSvc.h"
#include <list>


//forward declarations
template <class TYPE> class SvcFactory;
class IFlux;  // interface
class FluxMgr;  // actual manager
class IParticlePropertySvc; 

//!  Service that implements the IFluxSvc interface, to return an IFlux object.
//!  FluxSvc handles the creation and interfacing with Flux objects.  
class FluxSvc : virtual public Service, virtual public IFluxSvc
{  
public:
    
    /// return pointer to a flux object
    StatusCode source(std::string name, IFlux*&);
    
    /// return a list of possible names
    std::list<std::string> fluxNames()const;
    
    /// add a new SpectrumFactory
    virtual void addFactory(std::string name, const ISpectrumFactory* factory );
    
    /// access to the local random engine 
    virtual HepRandomEngine* getEngine();

    /// pass a specific amount of time
    virtual void pass ( double t);

    /// create a set of display windows using rootplot.
    void rootDisplay(std::vector<char*> arguments);

    ///return the pointer to the current IFlux object
    IFlux* currentFlux();

    
    //------------------------------------------------------------------
    //  stuff required by a Service
    
    /// perform initializations for this service. 
    virtual StatusCode initialize ();
    
    /// perform the finalization, as required for a service.
    virtual StatusCode finalize ();
    
    /// Query interface
    virtual StatusCode queryInterface( const IID& riid, void** ppvUnknown );
    
protected: 
    
    /// Standard Constructor
    FluxSvc ( const std::string& name, ISvcLocator* al );
    
    /// destructor
    virtual ~FluxSvc ();
    
private:

    IParticlePropertySvc* m_partSvc;

    /// Allow SvcFactory to instantiate the service.
    friend class SvcFactory<FluxSvc>;
    
    FluxMgr * m_fluxMgr;
    /// the user-defined list of acceptable XML sources (from JobOptions.txt)
    //std::string m_source_lib;
    std::vector<std::string> m_source_lib;
    /// the default XML file name (from JobOptions.txt)
    std::string m_source_lib_default;
    /// set dtd to use.
    std::string m_dtd_file;
    /// the "current" flux object
    IFlux* m_currentFlux;
};


#endif // _H_FluxSvc
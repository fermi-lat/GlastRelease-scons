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

/*! Service that implements the IFluxSvc interface, to return an IFlux object  
*/
class FluxSvc : virtual public Service, virtual public IFluxSvc
{  
public:

    /// return pointer to a flux object
    StatusCode source(std::string name, IFlux*&);

    /// return a list of possible names
    std::list<std::string> fluxNames()const;

    /// add a new source
    virtual void addFactory(std::string name, const ISpectrumFactory* factory );


    //------------------------------------------------------------------
    //  stuff required by a Service
    
    /// perform initializations for this service. 
    virtual StatusCode initialize ();
    
    /// 
    virtual StatusCode finalize ();

      
   /// Query interface
   virtual StatusCode queryInterface( const IID& riid, void** ppvUnknown );

protected: 

    /// Standard Constructor
    FluxSvc ( const std::string& name, ISvcLocator* al );
    
    /// destructor
    virtual ~FluxSvc ();

private:
    // Allow SvcFactory to instantiate the service.
    friend class SvcFactory<FluxSvc>;

    FluxMgr * m_fluxMgr;
};


#endif // _H_FluxSvc
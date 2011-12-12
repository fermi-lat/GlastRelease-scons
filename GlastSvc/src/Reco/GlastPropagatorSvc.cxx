

/** 
 * @class GlastPropagatorSvc
 *
 * @brief This service provides an interface for users to obtain the common
 *        propagator tool. Several implementations of a propagator exist, this
 *        tool will determine which to use and provide a mechanism for access to it.
 * 
 * @author Tracy Usher
 *
 * $Header$
 */

#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include <iostream>
#include <algorithm>

#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/Reco/IPropagatorTool.h"

class GlastPropagatorSvc : public Service, virtual public IPropagatorSvc
{
public:

    GlastPropagatorSvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~GlastPropagatorSvc() {}
    
    StatusCode initialize();
    StatusCode finalize();

     // This returns a pointer to the singleton propagator object
    virtual IKalmanParticle* getPropagator() const {return m_PropagatorTool->getPropagator();}
 
    /// Provide service interface
    StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);

    /// return the service type
    const InterfaceID& type() const;
    
private:

    /// Type of fit to perform
    std::string      m_PropagatorType;

    /// Pointer to the propagator guts
    IPropagatorTool* m_PropagatorTool;
};




//static const SvcFactory<GlastPropagatorSvc> s_factory;
//const ISvcFactory& GlastPropagatorSvcFactory = s_factory;
DECLARE_SERVICE_FACTORY( GlastPropagatorSvc );


//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.
/// The property declared is used to determine which propagator to use

GlastPropagatorSvc::GlastPropagatorSvc(const std::string& name, 
                               ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{   
    declareProperty("PropagatorType", m_PropagatorType="G4Propagator");

    m_PropagatorTool = 0;

    return; 
}

StatusCode GlastPropagatorSvc::initialize()
{
    // Purpose: Initialize the service and the specific propagator tool
    // Inputs:  none
    // Output:  none
    
    StatusCode sc = StatusCode::SUCCESS;

    // Take care of stuff to initialize the service
    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());

    // Instantiate the concrete tool we want for the propagator
    IToolSvc* toolSvc = 0;
    if (sc = service("ToolSvc",toolSvc, true).isSuccess() )
    {
        // Which propagator to use?
        if (m_PropagatorType == "G4Propagator")
        {
            // Look for the G4PropagatorSvc service
            sc = toolSvc->retrieveTool("G4PropagatorTool", m_PropagatorTool);
            log << MSG::INFO << "Using Geant4 Particle Propagator" << endreq;
        }
        else
        {
            // Look for GismoGenerator Service
            sc = toolSvc->retrieveTool("RecoTool", m_PropagatorTool);
            log << MSG::INFO << "Using Gismo Particle Propagator" << endreq;
        }
    } else { 
        log << MSG::INFO << "ToolSvc not found" << endreq;
        return sc; 
    }   
    return sc;
}

StatusCode GlastPropagatorSvc::finalize()
{
    return StatusCode::SUCCESS;
}

StatusCode  GlastPropagatorSvc::queryInterface (const InterfaceID& riid, void **ppvIF)
{
    if (IPropagatorSvc::interfaceID() == riid) {
        *ppvIF = dynamic_cast<IPropagatorSvc*> (this);
        return StatusCode::SUCCESS;
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
}


const InterfaceID&  GlastPropagatorSvc::type () const {
    return IID_IPropagatorSvc;
}



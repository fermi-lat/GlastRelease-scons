
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "HepRepCorbaSvc.h"

DECLARE_SERVICE_FACTORY(HepRepCorbaSvc);

//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

HepRepCorbaSvc::HepRepCorbaSvc(const std::string& name, 
                               ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{   

    return; 
}

StatusCode HepRepCorbaSvc::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());

    IToolSvc* toolSvc = 0;
    if (sc = service("ToolSvc",toolSvc, true).isSuccess() )
    {
        sc = toolSvc->retrieveTool("RegisterCorba", m_xmlTool);
        if (sc.isSuccess()) {
            log << MSG::INFO << "Retrieved DetDisplay" << endreq;
        } else {
            log << MSG::ERROR << "Couldn't retrieve DetDisplay" << endreq;
        }

    } else { 
        log << MSG::INFO << "ToolSvc not found" << endreq;
        return sc; 
    } 

    log << MSG::INFO << "HepRepCorbaSvc Initialized" << endreq;
    return StatusCode::SUCCESS;

}

StatusCode HepRepCorbaSvc::finalize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "HepRepCorbaSvc finalize called" << endreq;
    return StatusCode::SUCCESS;
}

// queryInterface

//StatusCode HepRepCorbaSvc::queryInterface (const InterfaceID& riid, void **ppvIF)
//{
//    if (IID_IHepRepCorbaSvc == riid) {
//        *ppvIF = dynamic_cast<IHepRepCorbaSvc*> (this);
//        return StatusCode::SUCCESS;
//    }
//    else {
//        return Service::queryInterface (riid, ppvIF);
//    }
//}

// access the type of this service

//const InterfaceID&  HepRepCorbaSvc::type () const {
//    return IID_IHepRepCorbaSvc;
//}


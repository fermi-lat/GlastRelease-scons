
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "DetDisplaySvc.h"

DECLARE_SERVICE_FACTORY(DetDisplaySvc);

//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

DetDisplaySvc::DetDisplaySvc(const std::string& name, 
                               ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{   

    return; 
}

StatusCode DetDisplaySvc::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());

    IToolSvc* toolSvc = 0;
    if (sc = service("ToolSvc",toolSvc, true).isSuccess() )
    {
        sc = toolSvc->retrieveTool("DetectorDisplay", m_detDisplayTool);
        if (sc.isSuccess()) {
            log << MSG::INFO << "Retrieved DetectorDisplay" << endreq;
        } else {
            log << MSG::ERROR << "Couldn't retrieve DetDisplay" << endreq;
        }
        sc = toolSvc->retrieveTool("MCdisplay", m_mcTool);
        if (sc.isSuccess()) {
            log << MSG::INFO << "Retrieved MCdisplay" << endreq;
        } else {
            log << MSG::ERROR << "Couldn't retrieve MCdisplay" << endreq;
        }
        sc = toolSvc->retrieveTool("StripDisplay", m_stripTool);
        if (sc.isSuccess()) {
            log << MSG::INFO << "Retrieved StripDisplay" << endreq;
        } else {
            log << MSG::ERROR << "Couldn't retrieve StripDisplay" << endreq;
        }

    } else { 
        log << MSG::INFO << "ToolSvc not found" << endreq;
        return sc; 
    } 

    log << MSG::INFO << "DetDisplaySvc Initialized" << endreq;
    return StatusCode::SUCCESS;

}

StatusCode DetDisplaySvc::finalize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "DetDisplaySvc finalize called" << endreq;
    return StatusCode::SUCCESS;
}

// queryInterface

//StatusCode  DetDisplaySvc::queryInterface (const InterfaceID& riid, void **ppvIF)
//{
 //   if (IID_IDetDisplaySvc == riid) {
//        *ppvIF = dynamic_cast<IDetDisplaySvc*> (this);
//        return StatusCode::SUCCESS;
//    }
//    else {
//        return Service::queryInterface (riid, ppvIF);
//    }
//}

// access the type of this service

//const InterfaceID&  DetDisplaySvc::type () const {
//    return IID_IDetDisplaySvc;
//}


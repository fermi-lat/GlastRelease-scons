
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "HepRepXmlSvc.h"

DECLARE_SERVICE_FACTORY(HepRepXmlSvc);

//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

HepRepXmlSvc::HepRepXmlSvc(const std::string& name, 
                               ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{   

    return; 
}

StatusCode HepRepXmlSvc::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());

    IToolSvc* toolSvc = 0;
    if (sc = service("ToolSvc",toolSvc, true).isSuccess() )
    {
        sc = toolSvc->retrieveTool("RegisterXml", m_xmlTool);
        if (sc.isSuccess()) {
            log << MSG::INFO << "Retrieved DetDisplay" << endreq;
        } else {
            log << MSG::ERROR << "Couldn't retrieve DetDisplay" << endreq;
        }

    } else { 
        log << MSG::INFO << "ToolSvc not found" << endreq;
        return sc; 
    } 

    log << MSG::INFO << "HepRepXmlSvc Initialized" << endreq;
    return StatusCode::SUCCESS;

}

StatusCode HepRepXmlSvc::finalize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "HepRepXmlSvc finalize called" << endreq;
    return StatusCode::SUCCESS;
}

// queryInterface

//StatusCode HepRepXmlSvc::queryInterface (const InterfaceID& riid, void **ppvIF)
//{
//    if (IID_IHepRepXmlSvc == riid) {
//        *ppvIF = dynamic_cast<IHepRepXmlSvc*> (this);
//        return StatusCode::SUCCESS;
//    }
//    else {
//        return Service::queryInterface (riid, ppvIF);
//    }
//}

// access the type of this service

//const InterfaceID&  HepRepXmlSvc::type () const {
//    return IID_IHepRepXmlSvc;
//}


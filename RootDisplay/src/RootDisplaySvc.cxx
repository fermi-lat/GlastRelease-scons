
#include "RootDisplaySvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IToolSvc.h"

DECLARE_SERVICE_FACTORY(RootDisplaySvc);

//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

RootDisplaySvc::RootDisplaySvc(const std::string& name, 
                               ISvcLocator* pSvcLocator) :
Service(name, pSvcLocator)
{   

    return; 
}

StatusCode RootDisplaySvc::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());

    IToolSvc* toolSvc = 0;
    if (sc = service("ToolSvc",toolSvc, true).isSuccess() )
    {
        sc = toolSvc->retrieveTool("RootMenu", m_rootMenuTool);
        if (sc.isSuccess()) {
            log << MSG::INFO << "Retrieved RootMenu" << endreq;
        } else {
            log << MSG::ERROR << "Couldn't retrieve RootMenu" << endreq;
        }
    } else { 
        log << MSG::INFO << "ToolSvc not found" << endreq;
        return sc; 
    } 

    log << MSG::INFO << "RootDisplaySvc Initialized" << endreq;
    return StatusCode::SUCCESS;

}

StatusCode RootDisplaySvc::finalize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "RootDisplaySvc finalize called" << endreq;
    return StatusCode::SUCCESS;
}

// queryInterface

//StatusCode  RootDisplaySvc::queryInterface (const InterfaceID& riid, void **ppvIF)
//{
  //  if (IID_IRootDisplaySvc == riid) {
  //      *ppvIF = dynamic_cast<IRootDisplaySvc*> (this);
  //      return StatusCode::SUCCESS;
  //  }
   // else {
  //      return Service::queryInterface (riid, ppvIF);
   // }
//}

// access the type of this service

//const InterfaceID&  RootDisplaySvc::type () const {
//    return IID_IRootDisplaySvc;
//}


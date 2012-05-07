
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "HepRepSvc/IHepRepSvc.h"

#include "HepRepXmlSvc.h"
#include "XmlStreamer.h"

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

    // Make the connection between corba and HepRepSvc here
    IHepRepSvc* hepRepSvc = 0;
    if (sc = service("HepRepSvc", hepRepSvc, true).isSuccess())
    {
        log << MSG::INFO << "Register corba server ..." << endreq;
  
        XmlStreamer* xml = new XmlStreamer();
        hepRepSvc->addStreamer("xml",xml);
    }
    else
    {
        log <<  MSG::INFO << "Could not locate HepRepSvc" << endreq;
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


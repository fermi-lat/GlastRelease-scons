// File and Version Information:
//      $Header$
//
// Description:
//      EventCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define EVENTCNV_CPP 

#include "GaudiKernel/RegistryEntry.h"
#include "EventCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/TopLevel/Event.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<EventCnv> s_factory;
const ICnvFactory& EventCnvFactory = s_factory;

EventCnv::EventCnv(ISvcLocator* svc)
: BaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event", objType(), "PASS");
}


StatusCode EventCnv::createObj(IOpaqueAddress* pAddress, 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an empty EventHeader on
    //   the TDS.
    refpObject = new Event::EventHeader();
    return StatusCode::SUCCESS;
};

StatusCode EventCnv::updateObj(int* iarray, Event::EventHeader* pEvent) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.
    MsgStream log(msgSvc(), "EventCnv");
    log << MSG::DEBUG << "EventCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}



// File and Version Information:
//      $Header$
//
// Description:
//      EventCnv is the concrete converter for the event header.
//
// Author(s):

#define EVENTCNV_CPP 

#include "GaudiKernel/RegistryEntry.h"
#include "EventCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "GlastEvent/TopLevel/Event.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<EventCnv> s_factory;
const ICnvFactory& EventCnvFactory = s_factory;


StatusCode EventCnv::updateObj(int* iarray, Event* pEvent)             {
  MsgStream log(msgSvc(), "EventCnv");
  log << MSG::DEBUG << "EventCnv::updateObj" << endreq;
  return StatusCode::SUCCESS;
}


EventCnv::EventCnv(ISvcLocator* svc)
: BaseCnv(classID(), svc)
//  : ItemCnv<Event>(svc)                                  
{
      declareObject("/Event", objType(), "PASS");
}



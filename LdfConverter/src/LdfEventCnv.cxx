// File and Version Information:
//      $Header$
//
// Description:
//      LdfEventCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfEventCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "LdfEventCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/TopLevel/Event.h"
#include "ldfReader/data/LatData.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<LdfEventCnv> s_factory;
const ICnvFactory& LdfEventCnvFactory = s_factory;

LdfEventCnv::LdfEventCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event", objType(), "PASS");
}


StatusCode LdfEventCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an empty EventHeader on
    //   the TDS.
    Event::EventHeader *header = new Event::EventHeader();
    refpObject = header;
    // Retrieve the LAT data for this event 
    unsigned int summary = ldfReader::LatData::instance()->summaryData().summary();
    // According to LAT Inter-module Communications Reference Manual
    // the event number is a 17 bit value, the 2 least significant are the
    // tag and the 15 most significant are in the eventNumber field of the summary
    unsigned eventNumber = EventSummary::eventNumber(summary);
    unsigned tag = EventSummary::tag(summary);
    //long actualEventNum = (tag << 15) | eventNumber;
    long actualEventNum = (eventNumber << 2) | tag;
    header->setEvent(actualEventNum);
    header->setRun(ldfReader::LatData::instance()->runId());

    return StatusCode::SUCCESS;
}

StatusCode LdfEventCnv::updateObj(int* , Event::EventHeader*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "LdfEventCnv");
    log << MSG::DEBUG << "LdfEventCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}


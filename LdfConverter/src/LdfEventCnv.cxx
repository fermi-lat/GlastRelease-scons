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
#include "facilities/Timestamp.h"
#include "astro/JulianDate.h"
#include "Event/Utilities/TimeStamp.h"

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
	
    // Also set the time in the Event::EventHeader
    facilities::Timestamp facTimeStamp(
             ldfReader::LatData::instance()->summaryData().timeSec(), 
             ldfReader::LatData::instance()->summaryData().timeNanoSec());
    double julTimeStamp = facTimeStamp.getJulian();
    astro::JulianDate julDate(julTimeStamp);
    // Find number of seconds since missionStart
    double tdsTimeInSeconds = julDate.seconds() - astro::JulianDate::missionStart().seconds();
    TimeStamp tdsTimeStamp(tdsTimeInSeconds);
    header->setTime(tdsTimeStamp);

    // Also set the LiveTime using the LDF.  LiveTime in LDF is in counts
    // where each count corresponds to 50 nsec.
    ldfReader::GemData ldfGem = ldfReader::LatData::instance()->getGem();
    // Convert to seconds.. where nanoseconds correspond to 1 billion seconds
    double ltSec = ldfGem.liveTime() * 50./1000000000.0;
    header->setLivetime(ltSec);


    // initialize the trigger to zero to promote processing if TriggerAlg is
    // run
    header->setTrigger(-1);

    return StatusCode::SUCCESS;
}

StatusCode LdfEventCnv::updateObj(int* , Event::EventHeader*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "LdfEventCnv");
    log << MSG::DEBUG << "LdfEventCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}


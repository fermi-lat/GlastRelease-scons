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

#include "ldfReader/data/LatData.h"
#include "facilities/Timestamp.h"
#include "astro/JulianDate.h"
#include "Event/Utilities/TimeStamp.h"
#include <math.h>

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
    MsgStream log(msgSvc(), "LdfEventCnv");
    Event::EventHeader *header = new Event::EventHeader();
    refpObject = header;
    // Retrieve the LAT data for this event 
    unsigned int summary = ldfReader::LatData::instance()->summaryData().summary();
    // According to LAT Inter-module Communications Reference Manual
    // the event number is a 17 bit value, the 2 least significant are the
    // tag and the 15 most significant are in the eventNumber field of the summary
    header->setEvent(ldfReader::LatData::instance()->getOsw().evtSequence());
    header->setRun(ldfReader::LatData::instance()->runId());
	
    // Also set the time in the Event::EventHeader
//    unsigned upperPpcTime = ldfReader::LatData::instance()->summaryData().upperPpcTimeBaseWord();
//    unsigned lowerPpcTime = ldfReader::LatData::instance()->summaryData().lowerPpcTimeBaseWord();
//    const double sixteenMHz = 16000000.;
    // To handle th 64-bit value - we separate the computation
    // The upper 32 bits would have to be shifted by 31 (or multiplied)
    // by 2^32 - we divide this by 16000000 to get the upperMultiplier
//    const double upperMultiplier = 268.435456; 
//    double lower = lowerPpcTime / sixteenMHz;
//    double upper = upperPpcTime * upperMultiplier;
//    double ppcSeconds =  upper + lower;
//    float ppcWholeSeconds = floor(ppcSeconds);
//    double frac = ppcSeconds - ppcWholeSeconds;
//    int ppcNanoSec = frac / 0.000000001;
//    log << MSG::DEBUG << "ppcSeconds: " << ppcSeconds << " frace = " << frac << endreq;
 //   log << MSG::DEBUG << "ppcNanoSec: " << ppcNanoSec << endreq;

    // To eliminate duplicate times, we make use of the PPC time as suggested 
    // by Jim Panetta.  Here we ignore the fractional seconds from the "real"
    // absolute time and add in the fractional seconds as computed by the PPC
    // word.  This should produce a unique time stamp for each event.
//    log << MSG::DEBUG << "seconds: " << ldfReader::LatData::instance()->summaryData().timeSec() << endreq;
//    log << MSG::DEBUG << "nanosec: " << ldfReader::LatData::instance()->summaryData().timeNanoSec() << endreq;
//    facilities::Timestamp facTimeStamp(
//             ldfReader::LatData::instance()->summaryData().timeSec(), 
//             ppcNanoSec);
//    double julTimeStamp = facTimeStamp.getJulian();
//    astro::JulianDate julDate(julTimeStamp);
    // Find number of seconds since missionStart
//    double tdsTimeInSeconds = julDate.seconds() - astro::JulianDate::missionStart().seconds();
//    TimeStamp tdsTimeStamp(tdsTimeInSeconds);
    TimeStamp tdsTimeStamp(ldfReader::LatData::instance()->getTimeInSecTds());
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


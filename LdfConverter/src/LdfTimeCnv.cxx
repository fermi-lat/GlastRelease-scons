// File and Version Information:
//      $Header$
//
// Description:
//      LdfTimeCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfTimeCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "LdfTimeCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "ldfReader/data/LatData.h"
#include "LdfEvent/LdfTime.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<LdfTimeCnv> s_factory;
const ICnvFactory& LdfTimeCnvFactory = s_factory;

LdfTimeCnv::LdfTimeCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event/Time", objType(), "PASS");
}


StatusCode LdfTimeCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will store Ebf time on the TDS.
    LdfEvent::LdfTime *time = new LdfEvent::LdfTime();
    refpObject = time;
    // Retrieve the LAT data for this event 
    unsigned int summary = ldfReader::LatData::instance()->summaryData().summary();
    time->initialize(ldfReader::LatData::instance()->summaryData().timeSec(), 
        ldfReader::LatData::instance()->summaryData().timeNanoSec(),
        ldfReader::LatData::instance()->summaryData().upperPpcTimeBaseWord(),
        ldfReader::LatData::instance()->summaryData().lowerPpcTimeBaseWord());
    return StatusCode::SUCCESS;
}

StatusCode LdfTimeCnv::updateObj(int* , LdfEvent::LdfTime*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "LdfTimeCnv");
    log << MSG::DEBUG << "LdfTimeCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}


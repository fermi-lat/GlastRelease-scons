// File and Version Information:
//      $Header$
//
// Description:
//      LdfEventSummaryCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfEventSummaryCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "LdfEventSummaryCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/TopLevel/Event.h"
#include "ldfReader/data/LatData.h"
#include "LdfEvent/EventSummaryData.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<LdfEventSummaryCnv> s_factory;
const ICnvFactory& LdfEventSummaryCnvFactory = s_factory;

LdfEventSummaryCnv::LdfEventSummaryCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event/EventSummary", objType(), "PASS");
}


StatusCode LdfEventSummaryCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an empty EventHeader on
    //   the TDS.
    // Retrieve the LAT data for this event 
    unsigned int ebfSummary = ldfReader::LatData::instance()->summaryData().summary();
    LdfEvent::EventSummaryData *summary = new LdfEvent::EventSummaryData(ebfSummary);
    summary->initEventFlags(ldfReader::LatData::instance()->getEventFlags());


    // Retrieve the contribution lengths and store them in the EventSummaryData
    unsigned long gemLen = ldfReader::LatData::instance()->getGem().lenInBytes();
    unsigned long oswLen = ldfReader::LatData::instance()->getOsw().lenInBytes();
    unsigned long errLen = ldfReader::LatData::instance()->getErr().lenInBytes();
    unsigned long diagLen = ldfReader::LatData::instance()->diagnostic()->lenInBytes();
    unsigned long aemLen = ldfReader::LatData::instance()->getAem().lenInBytes();
    unsigned long temLen[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    unsigned int i;
    for (i = 0; i < 16; i++) {
        if(ldfReader::LatData::instance()->getTower(i)) {
            temLen[i] = ldfReader::LatData::instance()->getTower(i)->getTem().lenInBytes();
        }
    }

    summary->initContribLen(temLen, gemLen, oswLen, errLen, diagLen, aemLen);

    refpObject = summary;

    return StatusCode::SUCCESS;
}

StatusCode LdfEventSummaryCnv::updateObj(int* , LdfEvent::EventSummaryData*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "LdfEventSummaryCnv");
    log << MSG::DEBUG << "LdfEventSummaryCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}


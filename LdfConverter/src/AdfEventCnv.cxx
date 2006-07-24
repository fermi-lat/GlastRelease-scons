// File and Version Information:
//      $Header$
//
// Description:
//      AdfEventCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define AdfEventCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "AdfEventCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "ldfReader/data/LatData.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<AdfEventCnv> s_factory;
const ICnvFactory& AdfEventCnvFactory = s_factory;

AdfEventCnv::AdfEventCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event/AncillaryData/AdfEvent", objType(), "PASS");
}


StatusCode AdfEventCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an AdfEvent on
    //   the TDS.
    MsgStream log(msgSvc(), "AdfEventCnv");
    AncillaryData::AdfEvent *adf = new AncillaryData::AdfEvent();
    refpObject  = adf;
    // Retrieve the LAT data for this event 
    const ldfReader::AdfData& adfLdf = ldfReader::LatData::instance()->getAdf();

    return StatusCode::SUCCESS;
}

StatusCode AdfEventCnv::updateObj(int* , AncillaryData::AdfEvent*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "AdfEventCnv");
    log << MSG::DEBUG << "AdfEventCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}


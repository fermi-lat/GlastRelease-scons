// File and Version Information:
//      $Header$
//
// Description:
//      AncillaryEventCnv is the concrete converter for the digi event header 
//      on the TDS /Event/AncillaryEvent
//
// Author(s):

#define AncillaryEventCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "AncillaryEventCnv.h"
#include "GaudiKernel/MsgStream.h"


// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<AncillaryEventCnv> s_factory;
const ICnvFactory& AncillaryEventCnvFactory = s_factory;

AncillaryEventCnv::AncillaryEventCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event/AncillaryEvent", objType(), "PASS");
}


StatusCode AncillaryEventCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an empty EventHeader on
    //   the TDS.
    DataObject *digi = new DataObject();
    refpObject = digi;
    return StatusCode::SUCCESS;
};

StatusCode AncillaryEventCnv::updateObj(int* , DataObject* ) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "AncillaryEventCnv");
    log << MSG::DEBUG << "AncillaryEventCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}


// File and Version Information:
//      $Header$
//
// Description:
//      MetaEventCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define MetaEventCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "MetaEventCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "ldfReader/data/LatData.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<MetaEventCnv> s_factory;
const ICnvFactory& MetaEventCnvFactory = s_factory;

MetaEventCnv::MetaEventCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event/Gem path on the TDS.
    declareObject("/Event/MetaEvent", objType(), "PASS");
}


StatusCode MetaEventCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will store LDF Gem on the TDS.
    const lsfDataStore::MetaEvent& lsfMetaEvent = ldfReader::LatData::instance()->getMetaEvent();
    LsfEvent::MetaEvent *me = new LsfEvent::MetaEvent(lsfMetaEvent);
    refpObject = me;

    return StatusCode::SUCCESS;
}

StatusCode MetaEventCnv::updateObj(int* , LsfEvent::MetaEvent*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "MetaEventCnv");
    log << MSG::DEBUG << "MetaEventCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}


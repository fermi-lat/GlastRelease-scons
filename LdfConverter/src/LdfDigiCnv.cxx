// File and Version Information:
//      $Header$
//
// Description:
//      LdfDigiCnv is the concrete converter for the digi event header on the TDS /Event/Digi
//
// Author(s):

#define LdfDigiCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "LdfDigiCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "Event/TopLevel/DigiEvent.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<LdfDigiCnv> s_factory;
const ICnvFactory& LdfDigiCnvFactory = s_factory;

LdfDigiCnv::LdfDigiCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event/Digi", objType(), "PASS");
}


StatusCode LdfDigiCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an empty EventHeader on
    //   the TDS.
    Event::DigiEvent *digi = new Event::DigiEvent();
    // Set fromMc to false
    digi->initialize(false);
    refpObject = digi;
    return StatusCode::SUCCESS;
};

StatusCode LdfDigiCnv::updateObj(int* , Event::EventHeader* ) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "LdfDigiCnv");
    log << MSG::DEBUG << "LdfDigiCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}


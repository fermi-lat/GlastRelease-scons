// File and Version Information:
//      $Header$
//
// Description:
//      LdfErrorCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfErrorCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "LdfErrorCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "ldfReader/data/LatData.h"

#include "LdfEvent/ErrorData.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<LdfErrorCnv> s_factory;
const ICnvFactory& LdfErrorCnvFactory = s_factory;

LdfErrorCnv::LdfErrorCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event/Error", objType(), "PASS");
}


StatusCode LdfErrorCnv::createObj(IOpaqueAddress* , 
                                    DataObject*& refpObject) {
    // Purpose and Method:  This converter will store LDF diagnostic data on the TDS

    LdfEvent::ErrorData *errCol = new LdfEvent::ErrorData();
    refpObject = errCol; 
    
    // Retrieve the LAT data for this event and print out its contents
    ldfReader::LatData* myLatData = ldfReader::LatData::instance();
    unsigned int i;
    for (i = 0; i < 16; i++) {
        if(myLatData->getTower(i)) {
            if (!myLatData->getTower(i)->getTem().exist()) continue;

            const ldfReader::ErrData* errLdf = myLatData->getTower(i)->getTem().getErr();

            if (!errLdf->exist()) continue;
            unsigned short towerId = myLatData->getTower(i)->getTowerId();
            LdfEvent::TowerErrorData err(towerId, errLdf->cal(), errLdf->tkr(), errLdf->phs(), errLdf->tmo());
            errCol->addTowerError(err);
        }
    }
    return StatusCode::SUCCESS;
}

StatusCode LdfErrorCnv::updateObj(int* , LdfEvent::ErrorData* ) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "LdfErrorCnv");
    log << MSG::DEBUG << "LdfErrorCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}


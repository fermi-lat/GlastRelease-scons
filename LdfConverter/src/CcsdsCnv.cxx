// File and Version Information:
//      $Header$
//
// Description:
//      CcsdsCnv is the concrete converter for the CCSDS items on TDS 
//
// Author(s):

#define CcsdsCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "CcsdsCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "ldfReader/data/LatData.h"
#include "lsfData/LsfCcsds.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<CcsdsCnv> s_factory;
const ICnvFactory& CcsdsCnvFactory = s_factory;

CcsdsCnv::CcsdsCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event/Ccsds path on the TDS.
    declareObject("/Event/Ccsds", objType(), "PASS");
}


StatusCode CcsdsCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will store LDF Gem on the TDS.
    const lsfData::LsfCcsds& lsfCcsds = ldfReader::LatData::instance()->getCcsds();
    LsfEvent::LsfCcsds *me = new LsfEvent::LsfCcsds(lsfCcsds);
    refpObject = me;

    return StatusCode::SUCCESS;
}

StatusCode CcsdsCnv::updateObj(int* , LsfEvent::LsfCcsds*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "CcsdsCnv");
    log << MSG::DEBUG << "CcsdsCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}


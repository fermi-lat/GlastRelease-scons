// File and Version Information:
//      $Header$
//
// Description:
//      CcsdsCnv is the concrete converter for the CCSDS items on TDS 
//
// Author(s):

#define EbfCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "EbfCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "ldfReader/data/LatData.h"
#include "lsfData/Ebf.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<EbfCnv> s_factory;
const ICnvFactory& EbfCnvFactory = s_factory;

EbfCnv::EbfCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event/Filter/Ebf path on the TDS.
    declareObject("/Event/Filter/Ebf", objType(), "PASS");
}


StatusCode EbfCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will store LDF Gem on the TDS.
    MsgStream log(msgSvc(), "EbfCnv");
    log << MSG::DEBUG << "EbfCnv::createObj" << endreq;
    const lsfData::Ebf& lsfEbf = ldfReader::LatData::instance()->getEbf();
    EbfWriterTds::Ebf *me = new EbfWriterTds::Ebf();
    unsigned int len;
    char *data = lsfEbf.get(len);
    log << MSG::INFO << "EBF len " << len  << " first byte "
        << std::hex << (unsigned int)(data[0]) << std::dec << endreq;
    me->set(data,len);
    refpObject = me;

    return StatusCode::SUCCESS;
}

StatusCode EbfCnv::updateObj(int* , EbfWriterTds::Ebf*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "EbfCnv");
    log << MSG::DEBUG << "EbfCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}

// File and Version Information:
//      $Header$
//
// Description:
//      DigiEventCnv is the concrete converter for the digi event header.
//
// Author(s):


#define CNV_DIGIEVENTCNV_CPP 

#include "GaudiKernel/CnvFactory.h"
#include "DigiEventCnv.h"
#include "Event/TopLevel/DigiEvent.h"

static const char* rcsid = "$Id$";

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<DigiEventCnv> s_factory;
const ICnvFactory& DigiEventCnvFactory = s_factory;

DigiEventCnv::DigiEventCnv(ISvcLocator* svc)
: BaseCnv(classID(), svc)
{
  declareObject("/Event/Digi", objType(), "PASS");
}


StatusCode DigiEventCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject) {
    refpObject = new Event::DigiEvent();
    return StatusCode::SUCCESS;
};





#define CNV_DIGIEVENTCNV_CPP 

#include "GaudiKernel/CnvFactory.h"
#include "DigiEventCnv.h"
#include "GlastEvent/TopLevel/DigiEvent.h"

static const char* rcsid = "$Id$";

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<DigiEventCnv> s_factory;
const ICnvFactory& DigiEventCnvFactory = s_factory;


StatusCode DigiEventCnv::updateObj(int* iarray, DigiEvent* pEvent)          {
  return StatusCode::SUCCESS;
}

DigiEventCnv::DigiEventCnv(ISvcLocator* svc)
: BaseCnv(classID(), svc)
{
  declareObject("/Event/Digi", objType(), "PASS");
}



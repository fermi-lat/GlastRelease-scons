//------------------------------------------------------------------------------
//
// Implementation of class :  IrfEventCnv
////
//------------------------------------------------------------------------------
#define CNV_IrfEventCNV_CPP 

// Include files
#include "Gaudi/Kernel/CnvFactory.h"
#include "IrfEventCnv.h"
#include "GlastEvent/TopLevel/IrfEvent.h"

// RCS Id for identification of object version
static const char* rcsid = "$Id$";

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<IrfEventCnv> s_factory;
const ICnvFactory& IrfEventCnvFactory = s_factory;


StatusCode IrfEventCnv::updateObj(int* iarray, IrfEvent* pEvent)          {
  return StatusCode::SUCCESS;
}

/// Standard Constructor
IrfEventCnv::IrfEventCnv(ISvcLocator* svc)
: ItemCnv<IrfEvent>(svc)               
{
  declareObject("/Event/Irf", objType(), "PASS");
}


/// Standard Destructor
IrfEventCnv::~IrfEventCnv()   { 
}

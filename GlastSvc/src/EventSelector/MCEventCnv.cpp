// File and Version Information:
// $Header$
//
// Description:
// Concrete converter for the McEvent header on the TDS /Event/MC

#define CNV_MCEVENTCNV_CPP 

#include "GaudiKernel/CnvFactory.h"
#include "MCEventCnv.h"

// Event for creating the McEvent stuff
#include "Event/TopLevel/MCEvent.h"

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ObjectVector.h"

// RCS Id for identification of object version
static const char* rcsid = "$Id$";

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<MCEventCnv> s_factory;
const ICnvFactory& MCEventCnvFactory = s_factory;


MCEventCnv::MCEventCnv(ISvcLocator* svc)
: BaseCnv(classID(), svc)
{
  declareObject("/Event/MC", objType(), "PASS");
}


StatusCode MCEventCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)
{
    // Purpose and Method:  Converter just creates an empty instance of MCEvent
    //   for the TDS.  The data members will be initialized by other components.

    refpObject = new Event::MCEvent();
    StatusCode sc=StatusCode::SUCCESS;
    return sc;
}


const CLID& MCEventCnv::classID(){ return Event::MCEvent::classID();}

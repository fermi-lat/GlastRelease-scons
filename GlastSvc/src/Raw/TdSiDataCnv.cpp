// TdSiDataCnv.cpp: implementation of the TdSiDataCnv class.
//
//////////////////////////////////////////////////////////////////////

#define TdSiDataCnv_CPP

#include "Gaudi/Kernel/CnvFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/data/TdCsIData.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "src/EventSelector/IRFConverter.h"
#include "GlastSvc/Raw/TdSiDataCnv.h"

/*! Handles the conversion of TdSiData class currently uses the IRF
    information.
*/


//! Instantiation of a static factory class used by clients to create instances
static CnvFactory<TdSiDataCnv> s_factory;
const ICnvFactory& TdSiDataCnvFactory = s_factory;

//! Create the transient representation of an object.
StatusCode TdSiDataCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)   {

    IRFConverter* myConverter = new IRFConverter();
    
    m_detSvc->accept(*myConverter);
    
    // Make sure that the IRFConverter copy of the TdSiData is
    // not being deleted
    refpObject =  myConverter->getTdSiData();
    
  StatusCode status = updateObj(pAddress, refpObject);

  // Avoid the memory leak.
  delete myConverter;
  myConverter = 0;

  // Check to see if all is well with the IRFCOnverter copy of
  // TdSiData
  if (refpObject == NULL )   {
    delete refpObject;
    refpObject = 0;
    return StatusCode::FAILURE;
  }
  
  return StatusCode::SUCCESS;
}

//! Update object from scratch
StatusCode TdSiDataCnv::updateObj(IOpaqueAddress* pAddress, DataObject* pObject)   {
    TdSiData* allData;

    // Try to dynamic_cast pObject into a TdCsIData
    try {
        allData  = dynamic_cast<TdSiData*>(pObject);
    } catch(...) {
       return StatusCode::FAILURE;
    }

    // Read in the ACD data via the GlastDetector::accept method
    // This is now being dynamically declared so that we can keep
    // the reference for the TDS.
    IRFConverter* myConverter = new IRFConverter();
    
    m_detSvc->accept(*myConverter);


    allData =  myConverter->getTdSiData();

    delete myConverter;
    myConverter = 0;

    if (pObject == NULL )   {
        delete pObject;
        pObject = 0;
        return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
}

//! class ID
const CLID& TdSiDataCnv::classID()    {
  return CLID_TdSiDataCnv; 
}

//! Standard Constructor
TdSiDataCnv::TdSiDataCnv(ISvcLocator* svc) 
: BaseCnv(classID(), svc)
{

  declareObject("/Event/Raw/TdSiDatas", objType(), "PASS");
}

//! Standard Destructor
TdSiDataCnv::~TdSiDataCnv(){ }

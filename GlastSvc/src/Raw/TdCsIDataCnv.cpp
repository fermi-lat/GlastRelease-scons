// TdCsIDataCnv.cpp: implementation of the TdCsIDataCnv class.
//
//////////////////////////////////////////////////////////////////////

#define TdCsIDataCnv_CPP

#include "Gaudi/Kernel/CnvFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/data/TdCsIData.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "src/EventSelector/IRFConverter.h"
#include "GlastSvc/Raw/TdCsIDataCnv.h"

/*! Handles the conversion of TdCsIData class currently uses the IRF
    information.
*/


//! Instantiation of a static factory class used by clients to create instances
static CnvFactory<TdCsIDataCnv> s_factory;
const ICnvFactory& TdCsIDataCnvFactory = s_factory;

//! Create the transient representation of an object.
StatusCode TdCsIDataCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)   {

  IRFConverter* myConverter = new IRFConverter();
  m_detSvc->accept(*myConverter);

  refpObject = myConverter->getTdCsIData();

  if ( refpObject == NULL )   {
    StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

//! Update object from scratch
StatusCode TdCsIDataCnv::updateObj(IOpaqueAddress* pAddress, DataObject* pObject)   {

    // read in the ACD data via the GlastDetector::accept method
    // Done dynamically so we can keep it alive for the TDS.
    IRFConverter* myConverter = new IRFConverter();
    
    m_detSvc->accept(*myConverter);

    pObject = myConverter->getTdCsIData();

    // Avoid the memory leak
    if(myConverter) {
        delete myConverter;
        myConverter = 0;
    }

    if(pObject != NULL)
    {
        return StatusCode::SUCCESS;
    }
      
    return StatusCode::FAILURE;
}

//! class ID
const CLID& TdCsIDataCnv::classID()    {
  return CLID_TdCsIDataCnv; 
}

//! Standard Constructor
TdCsIDataCnv::TdCsIDataCnv(ISvcLocator* svc) 
: BaseCnv(classID(), svc)
{

  declareObject("/Event/Raw/TdCsIDatas", objType(), "PASS");
}

//! Standard Destructor
TdCsIDataCnv::~TdCsIDataCnv(){ }

// TdCsIDataCnv.cpp: implementation of the TdCsIDataCnv class.
//
//////////////////////////////////////////////////////////////////////

#define TdCsIDataCnv_CPP

#include "Gaudi/Kernel/CnvFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/Raw/TdCsIData.h"
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
  refpObject = new TdCsIData(9);//change to actual number
  StatusCode status = updateObj(pAddress, refpObject);
  if ( !status.isSuccess() )   {
    delete refpObject;
    refpObject = 0;
  }
  return status;
}

//! Update object from scratch
StatusCode TdCsIDataCnv::updateObj(IOpaqueAddress* pAddress, DataObject* pObject)   {
    TdCsIData* allData;

    // Try to dynamic_cast pObject into a TdCsIData
    try {
        allData  = dynamic_cast<TdCsIData*>(pObject);
    } catch(...) {
       return StatusCode::FAILURE;
    }

    // read in the ACD data via the GlastDetector::accept method
    IRFConverter myConverter;
      m_detSvc->accept(myConverter);

    // use the copyUp method to load all the crystals into the new object
    allData->copyUp(myConverter.getTdCsIData(), 9);

    if(allData != NULL)
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

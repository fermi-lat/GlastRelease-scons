// CsIDataCnv.cpp: implementation of the CsIDataCnv class.
//
//////////////////////////////////////////////////////////////////////

#define CsIDataCnv_CPP

#include "Gaudi/Kernel/CnvFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/Raw/CsIData.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "src/EventSelector/IRFConverter.h"
#include "GlastSvc/Raw/CsIDataCnv.h"

/*! Handles the conversion of CsIData class currently uses the IRF
    information.
*/


//! Instantiation of a static factory class used by clients to create instances
static CnvFactory<CsIDataCnv> s_factory;
const ICnvFactory& CsIDataCnvFactory = s_factory;

//! Create the transient representation of an object.
StatusCode CsIDataCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)   {
  refpObject = new CsIData(9);//change to actual number
  StatusCode status = updateObj(pAddress, refpObject);
  if ( !status.isSuccess() )   {
    delete refpObject;
    refpObject = 0;
  }
  return status;
}

//! Update object from scratch
StatusCode CsIDataCnv::updateObj(IOpaqueAddress* pAddress, DataObject* pObject)   {
    CsIData* allData;

    // Try to dynamic_cast pObject into a CsIData
    try {
        allData  = dynamic_cast<CsIData*>(pObject);
    } catch(...) {
       return StatusCode::FAILURE;
    }

    // read in the ACD data via the GlastDetector::accept method
    IRFConverter myConverter;
      m_detSvc->accept(myConverter);

    // use the copyUp method to load all the crystals into the new object
    allData->copyUp(myConverter.getCsIData(), 9);

    if(allData != NULL)
    {
        return StatusCode::SUCCESS;
    }
      
    return StatusCode::FAILURE;
}

//! class ID
const CLID& CsIDataCnv::classID()    {
  return CLID_CsIDataCnv; 
}

//! Standard Constructor
CsIDataCnv::CsIDataCnv(ISvcLocator* svc) 
: BaseCnv(classID(), svc)
{

  declareObject("/Event/Raw/CsIData", objType(), "PASS");
}

//! Standard Destructor
CsIDataCnv::~CsIDataCnv(){ }

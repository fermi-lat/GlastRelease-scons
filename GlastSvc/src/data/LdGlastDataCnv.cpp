// LdGlastDataCnv.cpp: implementation of the LdGlastDataCnv class.
//
//////////////////////////////////////////////////////////////////////

#define LdGlastDataCnv_CPP

#include "Gaudi/Kernel/CnvFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastSvc/../src/data/TdCsIData.h"
#include "GlastSvc/../src/data/TdSiData.h"
#include "GlastSvc/../src/data/TdVetoData.h"
#include "GlastSvc/../src/data/LdGlastData.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "src/EventSelector/IRFConverter.h"
#include "GlastSvc/data/LdGlastDataCnv.h"



//! Instantiation of a static factory class used by clients to create instances
static CnvFactory<LdGlastDataCnv> s_factory;
const ICnvFactory& LdGlastDataCnvFactory = s_factory;

// local static copy.
static IRFConverter myConverter;

//! Create the transient representation of an object.
StatusCode LdGlastDataCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)   {

    myConverter.clear();
    
    m_detSvc->accept(myConverter);
    
    // Make sure that the IRFConverter copy of the TdSiData is
    // not being deleted
    
    const SiData* si = myConverter.getSiData();
    const CsIData* csi = myConverter.getCsIData();
    const IVetoData* veto = myConverter.getVetoData();
    
    LdGlastData* glastData = new LdGlastData(csi,si,veto);
    
    refpObject =  glastData;
    
    
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
StatusCode LdGlastDataCnv::updateObj(IOpaqueAddress* pAddress, DataObject* pObject)   {
/*    TdSiData* allData;

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
    }*/
    return StatusCode::SUCCESS;
    
}

//! class ID
const CLID& LdGlastDataCnv::classID()    {
  return CLID_LdGlastDataCnv; 
}

//! Standard Constructor
LdGlastDataCnv::LdGlastDataCnv(ISvcLocator* svc) 
: BaseCnv(classID(), svc)
{

  declareObject("/Event/TdGlastData", objType(), "PASS");
}

//! Standard Destructor
LdGlastDataCnv::~LdGlastDataCnv(){ }


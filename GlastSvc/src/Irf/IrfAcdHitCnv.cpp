//------------------------------------------------------------------------------
//
// Implementation of class :  IrfAcdHitCnv
//
//
//------------------------------------------------------------------------------
#define IrfAcdHitCnv_CPP 

// Include files
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/Irf/IrfAcdHit.h"
#include "GaudiKernel/IGlastDetSvc.h"

#include "src/EventSelector/IRFConverter.h"

// Inherited sources
#include "GlastSvc/Irf/IrfAcdHitCnv.h"


// Instantiation of a static factory class used by clients to create instances
static CnvFactory<IrfAcdHitCnv> s_factory;
const ICnvFactory& IrfAcdHitCnvFactory = s_factory;

/// Create the transient representation of an object.
StatusCode IrfAcdHitCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)   {
  refpObject = new IrfAcdHitVector();
  StatusCode status = updateObj(pAddress, refpObject);
  if ( !status.isSuccess() )   {
    delete refpObject;
    refpObject = 0;
  }
  return status;
}

/// Update object from scratch
StatusCode IrfAcdHitCnv::updateObj(IOpaqueAddress* pAddress, DataObject* pObject)   {
  ObjectVector<IrfAcdHit>* tileList  = dynamic_cast<ObjectVector<IrfAcdHit>*>(pObject);
  if ( 0 != tileList )    {

      // read in the ACD data via the GlastDetector::accept method
      IRFConverter myConverter;
      m_detSvc->accept(myConverter);
      
      // iterate over the tiles and store in the pObject vector
      ObjectVector<IrfAcdHit>* acdData = myConverter.getIrfAcdHits();// Very important method.
      for (ObjectVector<IrfAcdHit>::const_iterator it = acdData->begin(); it != acdData->end(); it++) {
          IrfAcdHit * tile = new IrfAcdHit();
          tile->setEnergy((*it)->energy());
          tile->setId((*it)->id());
          tileList->push_back(tile);
      }

      return StatusCode::SUCCESS;
  }

  return StatusCode::FAILURE;
}

/// class ID
const CLID& IrfAcdHitCnv::classID()    {
  return ObjectVector<IrfAcdHit>::classID(); 
}

/// Standard Constructor
IrfAcdHitCnv::IrfAcdHitCnv(ISvcLocator* svc) 
: BaseCnv(classID(), svc)
{

  declareObject("/Event/Irf/IrfAcdHits", objType(), "PASS");
}

/// Standard Destructor
IrfAcdHitCnv::~IrfAcdHitCnv(){ }

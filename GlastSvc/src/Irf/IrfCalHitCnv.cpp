//------------------------------------------------------------------------------
//
// Implementation of class  :  IrfCalHitCnv
//
// Desciption               :  Gaudi converter for IrfCalHit
//------------------------------------------------------------------------------

#define IrfCalHitCnv_CPP 

// Include files
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/Irf/IrfCalHit.h"
#include "GaudiKernel/IGlastDetSvc.h"

#include "src/EventSelector/IRFConverter.h"

#include "GlastSvc/Irf/IrfCalHitCnv.h"

/*! Handles the conversion of IrfCalHit class currently uses the IRF
    information.
    
      Current Data Handled
        -minus end Response
        -plus end Response
        -Energy
*/


//! Instantiation of a static factory class used by clients to create instances
static CnvFactory<IrfCalHitCnv> s_factory;
const ICnvFactory& IrfCalHitCnvFactory = s_factory;

//! Create the transient representation of an object.
StatusCode IrfCalHitCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)   {
  refpObject = new IrfCalHitVector();
  StatusCode status = updateObj(pAddress, refpObject);
  if ( !status.isSuccess() )   {
    delete refpObject;
    refpObject = 0;
  }
  return status;
}

//! Update object from scratch
StatusCode IrfCalHitCnv::updateObj(IOpaqueAddress* pAddress, DataObject* pObject)   {
  ObjectVector<IrfCalHit>* calHitList  = dynamic_cast<ObjectVector<IrfCalHit>*>(pObject);
  if ( 0 != calHitList )    {

      // read in the ACD data via the GlastDetector::accept method
      IRFConverter myConverter;
      m_detSvc->accept(myConverter);
      
      // iterate over the tiles and store in the pObject vector
      ObjectVector<IrfCalHit>* mcCalData = myConverter.getIrfCalHits();// Change this
      for (ObjectVector<IrfCalHit>::const_iterator it = mcCalData->begin(); it != mcCalData->end(); it++) {
          IrfCalHit * mcCalHit = new IrfCalHit();
          mcCalHit->setEnergy((*it)->energy());
          mcCalHit->setMinusResponse((*it)->minusResponse());
          mcCalHit->setPlusResponse((*it)->plusResponse());
          calHitList->push_back(mcCalHit);
      }

      return StatusCode::SUCCESS;
  }

  return StatusCode::FAILURE;
}

//! class ID
const CLID& IrfCalHitCnv::classID()    {
  return ObjectVector<IrfCalHit>::classID(); 
}

//! Standard Constructor
IrfCalHitCnv::IrfCalHitCnv(ISvcLocator* svc) 
: BaseCnv(classID(), svc)
{

  declareObject("/Event/Irf/IrfCalHits", objType(), "PASS");
}

//! Standard Destructor
IrfCalHitCnv::~IrfCalHitCnv(){ }
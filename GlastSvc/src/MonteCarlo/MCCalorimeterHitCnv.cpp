// $Header$
//------------------------------------------------------------------------------
//
// Implementation of class  :  MCCalorimeterHitCnv
//
// Desciption               :  Gaudi converter for MCCalorimeterHit
//------------------------------------------------------------------------------

#define MCCalorimeterHitCnv_CPP 

// Include files
#include "Gaudi/Kernel/CnvFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/MonteCarlo/MCCalorimeterHit.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "src/EventSelector/IRFConverter.h"

#include "GlastSvc/MonteCarlo/MCCalorimeterHitCnv.h"

/*! Handles the conversion of MCCalorimeterHit class currently uses the IRF
    information.
    
      Current Data Handled
        -Left Response
        -Right Response
        -Energy
*/


//! Instantiation of a static factory class used by clients to create instances
static CnvFactory<MCCalorimeterHitCnv> s_factory;
const ICnvFactory& MCCalorimeterHitCnvFactory = s_factory;

//! Create the transient representation of an object.
StatusCode MCCalorimeterHitCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)   {
  refpObject = new MCCalorimeterHitVector();
  StatusCode status = updateObj(pAddress, refpObject);
  if ( !status.isSuccess() )   {
    delete refpObject;
    refpObject = 0;
  }
  return status;
}

//! Update object from scratch
StatusCode MCCalorimeterHitCnv::updateObj(IOpaqueAddress* pAddress, DataObject* pObject)   {
  ObjectVector<MCCalorimeterHit>* calHitList  = dynamic_cast<ObjectVector<MCCalorimeterHit>*>(pObject);
  if ( 0 != calHitList )    {

      // read in the ACD data via the GlastDetector::accept method
      IRFConverter myConverter;
      m_detSvc->accept(myConverter);
      
      // iterate over the tiles and store in the pObject vector
      ObjectVector<MCCalorimeterHit>* mcCalData = myConverter.getMCCalHits();// Change this
      for (ObjectVector<MCCalorimeterHit>::const_iterator it = mcCalData->begin(); it != mcCalData->end(); it++) {
          MCCalorimeterHit * mcCalHit = new MCCalorimeterHit();
          mcCalHit->setEnergy((*it)->energy());
          mcCalHit->setLeftResponse((*it)->leftResponse());
          mcCalHit->setRightResponse((*it)->rightResponse());
          calHitList->push_back(mcCalHit);
      }

      return StatusCode::SUCCESS;
  }

  return StatusCode::FAILURE;
}

//! class ID
const CLID& MCCalorimeterHitCnv::classID()    {
  return ObjectVector<MCCalorimeterHit>::classID(); 
}

//! Standard Constructor
MCCalorimeterHitCnv::MCCalorimeterHitCnv(ISvcLocator* svc) 
: BaseCnv(classID(), svc)
{

  declareObject("/Event/MC/MCCalorimeterHits", objType(), "PASS");
}

//! Standard Destructor
MCCalorimeterHitCnv::~MCCalorimeterHitCnv(){ }
// $Header$
//------------------------------------------------------------------------------
//
// Implementation of class :  ACDhitCnv
//
//
//------------------------------------------------------------------------------
#define MCACDhitCnv_CPP 

// Include files
#include "Gaudi/Kernel/CnvFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/MonteCarlo/MCACDHit.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "src/EventSelector/IRFConverter.h"

// Inherited sources
#include "GlastSvc/MonteCarlo/MCACDHitCnv.h"


// Instantiation of a static factory class used by clients to create instances
static CnvFactory<MCACDHitCnv> s_factory;
const ICnvFactory& MCACDHitCnvFactory = s_factory;

/// Create the transient representation of an object.
StatusCode MCACDHitCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)   {
  refpObject = new MCACDHitVector();
  StatusCode status = updateObj(pAddress, refpObject);
  if ( !status.isSuccess() )   {
    delete refpObject;
    refpObject = 0;
  }
  return status;
}

/// Update object from scratch
StatusCode MCACDHitCnv::updateObj(IOpaqueAddress* pAddress, DataObject* pObject)   {
  ObjectVector<MCACDHit>* tileList  = dynamic_cast<ObjectVector<MCACDHit>*>(pObject);
  if ( 0 != tileList )    {

      // read in the ACD data via the GlastDetector::accept method
      IRFConverter myConverter;
      m_detSvc->accept(myConverter);
      
      // iterate over the tiles and store in the pObject vector
      ObjectVector<MCACDHit>* acdData = myConverter.getACDTiles();// Very important method.
      for (ObjectVector<MCACDHit>::const_iterator it = acdData->begin(); it != acdData->end(); it++) {
          MCACDHit * tile = new MCACDHit();
          tile->setEnergy((*it)->energy());
          tile->setId((*it)->id());
          tileList->push_back(tile);
      }

      return StatusCode::SUCCESS;
  }

  return StatusCode::FAILURE;
}

/// class ID
const CLID& MCACDHitCnv::classID()    {
  return ObjectVector<MCACDHit>::classID(); 
}

/// Standard Constructor
MCACDHitCnv::MCACDHitCnv(ISvcLocator* svc) 
: BaseCnv(classID(), svc)
{

  //declareObject("/Event/ACDTile", objType(), "PASS");
  declareObject("/Event/MC/MCACDHits", objType(), "PASS");
}

/// Standard Destructor
MCACDHitCnv::~MCACDHitCnv(){ }

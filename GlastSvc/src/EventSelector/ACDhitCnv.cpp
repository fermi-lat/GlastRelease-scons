// $Header$
//------------------------------------------------------------------------------
//
// Implementation of class :  ACDhitCnv
//
//
//------------------------------------------------------------------------------
#define ACDhitCnv_CPP 

// Include files
#include "Gaudi/Kernel/CnvFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/Hits/ACDhit.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "IRFConverter.h"

// Inherited sources
#include "src/EventSelector/ACDhitCnv.h"


// Instantiation of a static factory class used by clients to create instances
static CnvFactory<ACDhitCnv> s_factory;
const ICnvFactory& ACDhitCnvFactory = s_factory;

/// Create the transient representation of an object.
StatusCode ACDhitCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)   {
  refpObject = new ACDhitVector();
  StatusCode status = updateObj(pAddress, refpObject);
  if ( !status.isSuccess() )   {
    delete refpObject;
    refpObject = 0;
  }
  return status;
}

/// Update object from scratch
StatusCode ACDhitCnv::updateObj(IOpaqueAddress* pAddress, DataObject* pObject)   {
  ObjectVector<ACDhit>* tileList  = dynamic_cast<ObjectVector<ACDhit>*>(pObject);
  if ( 0 != tileList )    {

      // read in the ACD data via the GlastDetector::accept method
      IRFConverter myConverter;
      m_detSvc->accept(myConverter);
      
      // iterate over the tiles and store in the pObject vector
      ObjectVector<ACDhit>* acdData = myConverter.getACDTiles();
      for (ObjectVector<ACDhit>::const_iterator it = acdData->begin(); it != acdData->end(); it++) {
          ACDhit * tile = new ACDhit();
          tile->setEnergy((*it)->energy());
          tile->setId((*it)->id());
          tileList->push_back(tile);
      }

      return StatusCode::SUCCESS;
  }

  return StatusCode::FAILURE;
}

/// class ID
const CLID& ACDhitCnv::classID()    {
  return ObjectVector<ACDhit>::classID(); 
}

/// Standard Constructor
ACDhitCnv::ACDhitCnv(ISvcLocator* svc) 
: BaseCnv(classID(), svc)
{

  //declareObject("/Event/ACDTile", objType(), "PASS");
  declareObject("/Event/MC/MCACDHits", objType(), "PASS");
}

/// Standard Destructor
ACDhitCnv::~ACDhitCnv(){ }

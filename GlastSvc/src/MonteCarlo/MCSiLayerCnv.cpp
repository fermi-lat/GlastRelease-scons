//------------------------------------------------------------------------------
//
// Implementation of class :  MCSiLayerCnv
//
//
//------------------------------------------------------------------------------
#define MCSiLayerCnv_CPP 

// Include files
#include "Gaudi/Kernel/CnvFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/MonteCarlo/MCSiLayer.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "src/EventSelector/IRFConverter.h"

// Inherited sources
#include "GlastSvc/MonteCarlo/MCSiLayerCnv.h"

/*! Converter that handles TRK data from the IRF file.
*/


//! Instantiation of a static factory class used by clients to create instances
static CnvFactory<MCSiLayerCnv> s_factory;
const ICnvFactory& MCSiLayerCnvFactory = s_factory;

//! Create the transient representation of an object.
StatusCode MCSiLayerCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)   {
  refpObject = new MCSiLayerVector();
  StatusCode status = updateObj(pAddress, refpObject);
  if ( !status.isSuccess() )   {
    delete refpObject;
    refpObject = 0;
  }
  return status;
}

/// Update object from scratch
StatusCode MCSiLayerCnv::updateObj(IOpaqueAddress* pAddress, DataObject* pObject)   {
  ObjectVector<MCSiLayer>* silayerList  = dynamic_cast<ObjectVector<MCSiLayer>*>(pObject);
  if ( 0 != silayerList )    {

      // read in the TKR data via the GlastDetector::accept method
      IRFConverter myConverter;
      m_detSvc->accept(myConverter);
      
      // iterate over the tiles and store in the pObject vector
      ObjectVector<MCSiLayer>* tkrData = myConverter.getMCTKRHits();// Very important method.
      for (ObjectVector<MCSiLayer>::const_iterator it = tkrData->begin(); it != tkrData->end(); it++) {
          MCSiLayer* layer = new MCSiLayer();
          layer->setMaxEnergy((*it)->MaxEnergy());
          layer->setId((*it)->id());
          // iterate over the hit strips
          for (ObjectVector<MCTKRHit>::const_iterator hit = ((*it)->getHits())->begin(); hit != ((*it)->getHits())->end(); hit++) {
              MCTKRHit* curHit = new MCTKRHit();
              curHit->setEnergy((*hit)->energy());
              curHit->setId((*hit)->id());
              curHit->setNoise((*hit)->noise());
              layer->addHit(curHit);
          }
          silayerList->push_back(layer);
      }

      return StatusCode::SUCCESS;
  }

  return StatusCode::FAILURE;
}

/// class ID
const CLID& MCSiLayerCnv::classID()    {
  return ObjectVector<MCSiLayer>::classID(); 
}

/// Standard Constructor
MCSiLayerCnv::MCSiLayerCnv(ISvcLocator* svc) 
: BaseCnv(classID(), svc)
{

  declareObject("/Event/MC/MCTKRHits", objType(), "PASS");
}

/// Standard Destructor
MCSiLayerCnv::~MCSiLayerCnv(){ }

//------------------------------------------------------------------------------
//
// Implementation of class :  IrfTkrLayerCnv
//
//
//------------------------------------------------------------------------------
#define IrfTkrLayerCnv_CPP 

// Include files
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/Irf/IrfTkrLayer.h"
#include "GaudiKernel/IGlastDetSvc.h"

#include "src/EventSelector/IRFConverter.h"

// Inherited sources
#include "GlastSvc/Irf/IrfTkrLayerCnv.h"

/*! Converter that handles TRK data from the IRF file.
*/


//! Instantiation of a static factory class used by clients to create instances
static CnvFactory<IrfTkrLayerCnv> s_factory;
const ICnvFactory& IrfTkrLayerCnvFactory = s_factory;

//! Create the transient representation of an object.
StatusCode IrfTkrLayerCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)   {
  refpObject = new IrfTkrLayerVector();
  StatusCode status = updateObj(pAddress, refpObject);
  if ( !status.isSuccess() )   {
    delete refpObject;
    refpObject = 0;
  }
  return status;
}

/// Update object from scratch
StatusCode IrfTkrLayerCnv::updateObj(IOpaqueAddress* pAddress, DataObject* pObject)   {
  ObjectVector<IrfTkrLayer>* silayerList  = dynamic_cast<ObjectVector<IrfTkrLayer>*>(pObject);
  if ( 0 != silayerList )    {

      // read in the TKR data via the GlastDetector::accept method
      IRFConverter myConverter;
      m_detSvc->accept(myConverter);
      
      // iterate over the tiles and store in the pObject vector
      ObjectVector<IrfTkrLayer>* tkrData = myConverter.getIrfTkrHits();// Very important method.
      for (ObjectVector<IrfTkrLayer>::const_iterator it = tkrData->begin(); it != tkrData->end(); it++) {
          IrfTkrLayer* layer = new IrfTkrLayer();
          layer->setMaxEnergy((*it)->MaxEnergy());
          layer->setId((*it)->id());
          // iterate over the hit strips
          for (ObjectVector<IrfTkrHit>::const_iterator hit = ((*it)->getHits())->begin(); hit != ((*it)->getHits())->end(); hit++) {
              IrfTkrHit* curHit = new IrfTkrHit();
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
const CLID& IrfTkrLayerCnv::classID()    {
  return ObjectVector<IrfTkrLayer>::classID(); 
}

/// Standard Constructor
IrfTkrLayerCnv::IrfTkrLayerCnv(ISvcLocator* svc) 
: BaseCnv(classID(), svc)
{
    declareObject("/Event/Irf/IrfTkrHits", objType(), "PASS");
}

/// Standard Destructor
IrfTkrLayerCnv::~IrfTkrLayerCnv(){ }

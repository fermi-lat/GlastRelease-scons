// File and Version Information:
// $Header$
//
// Description: Service for particle transport management
//
// Author(s):
//      T.Usher

#include "G4PropagatorSvc.h"
#include "G4ParticlePropagator.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "GaudiKernel/Service.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"

#include <cassert>

static SvcFactory<G4PropagatorSvc> g4prop_factory;

const ISvcFactory& G4PropagatorSvcFactory = g4prop_factory;

G4VUserDetectorConstruction* G4PropagatorSvc::UserDetector = 0;
IG4GeometrySvc::IdMap* G4PropagatorSvc::idmap = 0;
G4TransportationManager* G4PropagatorSvc::TransportationManager = 0;

G4PropagatorSvc::G4PropagatorSvc(const std::string& name, ISvcLocator* pSvcLocator) :
  Service(name, pSvcLocator)
{
  return;
}

G4PropagatorSvc::~G4PropagatorSvc()
{
}

StatusCode G4PropagatorSvc::initialize()
{
  // Purpose and Method:  Gaudi initialization routine. 
  // Inputs:  None
  // Outputs:  Gaudi StatusCode
  // Dependencies: None
  // Restrictions None 

  MsgStream log(msgSvc(), name());

  StatusCode sc = Service::initialize();
    
  sc = setProperties();

  // Recover the pointer to the GLAST Geant run manager. If it does not already
  // exist, then instantiate one.
  // We need this for using Geant to do the tracking
  // Get the Glast detector service 
  IGlastDetSvc* gdsv=0;
  if( service( "GlastDetSvc", gdsv).isFailure() ) 
    {
      log << MSG::ERROR << "Couldn't set up GlastDetSvc!" << endreq;
      return 0;
    }
  IG4GeometrySvc* gsv=0;
  if( service( "G4GeometrySvc", gsv).isFailure() ) 
    {
      log << MSG::ERROR << "Couldn't set up GlastDetSvc!" << endreq;
      return 0;
    }

   UserDetector = gsv->getDetector();
   TransportationManager = gsv->getTransportationManager();

   log << MSG::DEBUG << "G4 RunManager ready" << endreq;

  idmap = gsv->getIdMap();

  return sc;
}


StatusCode G4PropagatorSvc::finalize()
{
  // Purpose and Method:  Gaudi finalize routine. 
  // Inputs:  None
  // Outputs:  Gaudi StatusCode
  // Dependencies: None
  // Restrictions None 

  //If no run manager then delete the private version of the geometry
  //if (!RunManager::GetRunManager())
  //  {
      delete UserDetector;
  //  }

  StatusCode sc = Service::finalize();

  return sc;
}


// Query interface
StatusCode G4PropagatorSvc::queryInterface(const IID& riid, 
                                           void** ppvInterface)  
{
  // Purpose and Method:  Gaudi service query interface routine 
  // Inputs:  None
  // Outputs:  Gaudi StatusCode
  // Dependencies: None
  // Restrictions None 

  if ( IID_IPropagatorSvc.versionMatch(riid) ) 
    {
      *ppvInterface = (G4PropagatorSvc*)this;
    }
  else  
    {
      return Service::queryInterface(riid, ppvInterface);
    }
  addRef();
  return SUCCESS;
}

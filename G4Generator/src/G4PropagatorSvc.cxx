// File and Version Information:
// $Header$
//
// Description: Service for particle transport management
//
// Author(s):
//      T.Usher

#include "G4PropagatorSvc.h"

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/Property.h"

#include "G4TransportationManager.hh"

#include <cassert>


static SvcFactory<G4PropagatorSvc> g4_factory;

const ISvcFactory& G4PropagatorSvcFactory = g4_factory;

DetectorConstruction* G4PropagatorSvc::UserDetector = 0;
DetectorConstruction::IdMap* G4PropagatorSvc::idmap = 0;

//G4PropagatorSvc::G4PropagatorSvc(const std::string& name, ISvcLocator*
//pSvcLocator) : DataSvc(name, pSvcLocator)
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
  if (RunManager::GetRunManager())
    {
      const G4VUserDetectorConstruction* pUser = 
         RunManager::GetRunManager()->GetUserDetectorConstruction();
      const DetectorConstruction*        pDet = 
         dynamic_cast<const DetectorConstruction*>(pUser);
      UserDetector = const_cast<DetectorConstruction*>(pDet);
    }
  else
    {
      // Get the Glast detector service 
      IGlastDetSvc* gsv=0;
      if( service( "GlastDetSvc", gsv).isFailure() ) 
        {
          log << MSG::ERROR << "Couldn't set up GlastDetSvc!" << endreq;
          return 0;
        }

      // Retrieve the event data service
      IDataProviderSvc* eventSvc=0;
      if (service( "EventDataSvc", eventSvc, true ).isFailure())
        {
          log << MSG::ERROR << "Couldn't set up EventDataSvc!" << endreq;
          return 0;
        }

      UserDetector = new DetectorConstruction(gsv,eventSvc, "recon", std::cout);
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->SetWorldVolume(UserDetector->Construct());

      log << MSG::DEBUG << "G4 RunManager ready" << endreq;
    }

  idmap = UserDetector->idMap();

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
  if (!RunManager::GetRunManager())
    {
      delete UserDetector;
    }

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

// File and Version Information:
// $Header$
//
// Description: Service for particle transport management
//
// Author(s):
//      T.Usher

#include "G4Generator/IG4GeometrySvc.h"

#include "GaudiKernel/Service.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"

#include "G4VUserDetectorConstruction.hh"
#include "G4TransportationManager.hh"
#include "DetectorConstruction.h"

#include <cassert>

/** 
 * @class G4GeometrySvc
 *
 * @brief A Gaudi Service for setting up the geometry for Geant4 and the Propagator
 *
 * @author Tracy Usher
 *
 */
class G4GeometrySvc : public Service, virtual public IG4GeometrySvc
{
 public: 

  G4GeometrySvc(const std::string& name, ISvcLocator* pSvcLocator);

  virtual ~G4GeometrySvc();

  virtual StatusCode initialize();

  virtual StatusCode finalize();
        
  /// queryInterface - for implementing a Service this is necessary
  StatusCode queryInterface(const IID& riid, void** ppvUnknown);

  /// Return pointer to the constructed Detector
  virtual G4VUserDetectorConstruction* getDetector();

  /// Return pointer to the G4 Transportation Manager
  virtual G4TransportationManager* getTransportationManager();

  /// Return pointer to the ID map
  virtual IG4GeometrySvc::IdMap* getIdMap();
  
 private:
  /// This is needed by propagator for decoding volume identifiers
  G4VUserDetectorConstruction* UserDetector;
  G4TransportationManager* TransportationManager;
  IG4GeometrySvc::IdMap* idmap;

  /// the geometry level of details
  std::string m_geometryMode;

  // Allow SvcFactory to instantiate the service.
  friend class SvcFactory<G4GeometrySvc>;
};


static SvcFactory<G4GeometrySvc> g4_factory;
const ISvcFactory& G4GeometrySvcFactory = g4_factory;

G4GeometrySvc::G4GeometrySvc(const std::string& name, ISvcLocator* pSvcLocator) :
  Service(name, pSvcLocator), UserDetector(0), TransportationManager(0), idmap(0)
{
  declareProperty("geometryMode", m_geometryMode="recon");

  return;
}

G4GeometrySvc::~G4GeometrySvc()
{
  //Clean up the geometry if it is there
  if (UserDetector)
  {
	  DetectorConstruction* UserDetectorConstruction = dynamic_cast<DetectorConstruction*>(UserDetector);
      delete UserDetectorConstruction;
  }
}

StatusCode G4GeometrySvc::initialize()
{
  // Purpose and Method:  Gaudi initialization routine. 
  // Inputs:  None
  // Outputs:  Gaudi StatusCode
  // Dependencies: None
  // Restrictions None 

  MsgStream log(msgSvc(), name());
  
  log << MSG::INFO << "G4GeometrySvc Initialization started" << endreq;
  
  StatusCode sc = Service::initialize();
    
  sc = setProperties();
    
  // Recover the pointer to the GLAST Geant run manager. If it does not already
  // exist, then instantiate one.
  // We need this for using Geant to do the tracking
  // Get the Glast detector service 
  IGlastDetSvc* gsv=0;
  if( service( "GlastDetSvc", gsv, true).isFailure() ) 
    {
      log << MSG::ERROR << "Couldn't set up GlastDetSvc!" << endreq;
      return StatusCode::FAILURE;
    }

  // Retrieve the event data service
  IDataProviderSvc* eventSvc=0;
  if (service( "EventDataSvc", eventSvc, true ).isFailure())
    {
      log << MSG::ERROR << "Couldn't set up EventDataSvc!" << endreq;
      return StatusCode::FAILURE;
    }

  DetectorConstruction* UserDetectorConstruction = new DetectorConstruction(gsv, eventSvc, m_geometryMode, std::cout);

  TransportationManager = G4TransportationManager::GetTransportationManager();

  TransportationManager->GetNavigatorForTracking()->SetWorldVolume(UserDetectorConstruction->Construct());

  UserDetector          = UserDetectorConstruction;
  idmap                 = UserDetectorConstruction->idMap();
  
  log << MSG::INFO << "G4GeometrySvc Initialization completed successfully" << endreq;

  return sc;
}


StatusCode G4GeometrySvc::finalize()
{
  // Purpose and Method:  Gaudi finalize routine. 
  // Inputs:  None
  // Outputs:  Gaudi StatusCode
  // Dependencies: None
  // Restrictions None 

  StatusCode sc = Service::finalize();

  return sc;
}

StatusCode G4GeometrySvc::queryInterface(const IID& riid, void** ppvInterface)  
{
  // Purpose and Method:  Gaudi service query interface routine 
  // Inputs:  None
  // Outputs:  Gaudi StatusCode
  // Dependencies: None
  // Restrictions None 

  if ( IID_IG4GeometrySvc.versionMatch(riid) ) 
    {
      *ppvInterface = dynamic_cast<IG4GeometrySvc*>(this);
    }
  else  
    {
      return Service::queryInterface(riid, ppvInterface);
    }
  return SUCCESS;
}

/// Return pointer to the constructed Detector
G4VUserDetectorConstruction* G4GeometrySvc::getDetector() 
{
	return UserDetector;
}

/// Return pointer to the G4 Transportation Manager
G4TransportationManager* G4GeometrySvc::getTransportationManager() 
{
	return TransportationManager;
}

/// Return pointer to the ID map
IG4GeometrySvc::IdMap* G4GeometrySvc::getIdMap() 
{
	return idmap;
}

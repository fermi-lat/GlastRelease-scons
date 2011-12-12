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

// needed to supress warning messages, same as in DetectorConstruction.h 
#ifdef WIN32  
# include <float.h>
#endif

#include "G4VUserDetectorConstruction.hh"
#include "G4TransportationManager.hh"

#include "DetectorConstruction.h"
#include "LocalMagneticFieldDes.h"

//#include <cassert>

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
  StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);

  /// Return pointer to the constructed Detector
  virtual G4VUserDetectorConstruction* getDetector() {return m_UserDetector;}

  /// Return pointer to the G4 Transportation Manager
  virtual G4TransportationManager* getTransportationManager() {return m_TransportationManager;}

  /// Return identifier associated with a given volume
  virtual idents::VolumeIdentifier getVolumeIdent(const G4VPhysicalVolume* volume);
  
 private:
  /// This is needed by propagator for decoding volume identifiers
  G4VUserDetectorConstruction* m_UserDetector;
  G4TransportationManager*     m_TransportationManager;
  DetectorConstruction::IdMap* m_idmap;

  /// the geometry level of details
  std::string m_geometryMode;

  LocalMagneticFieldDes m_magFieldDes;

  // Allow SvcFactory to instantiate the service.
  friend class SvcFactory<G4GeometrySvc>;
};


//static SvcFactory<G4GeometrySvc> g4_factory;
//const ISvcFactory& G4GeometrySvcFactory = g4_factory;
DECLARE_SERVICE_FACTORY(G4GeometrySvc);

G4GeometrySvc::G4GeometrySvc(const std::string& name, ISvcLocator* pSvcLocator) :
  Service(name, pSvcLocator), m_UserDetector(0), m_TransportationManager(0), m_idmap(0)
{
  declareProperty("geometryMode", m_geometryMode="recon");
  declareProperty("magneticFieldVolumn", m_magFieldDes.m_magFieldVol="");
  declareProperty("magneticFieldX", m_magFieldDes.m_magFieldX=0);
  declareProperty("magneticFieldY", m_magFieldDes.m_magFieldY=0);
  declareProperty("magneticFieldZ", m_magFieldDes.m_magFieldZ=0);
  return;
}

G4GeometrySvc::~G4GeometrySvc()
{
  //Clean up the geometry if it is there
  if (m_UserDetector)
  {
	  DetectorConstruction* UserDetectorConstruction = dynamic_cast<DetectorConstruction*>(m_UserDetector);
      //delete UserDetectorConstruction;
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

  DetectorConstruction* UserDetectorConstruction = new DetectorConstruction(gsv, eventSvc, m_geometryMode, std::cout, m_magFieldDes);

  m_TransportationManager = G4TransportationManager::GetTransportationManager();

  m_TransportationManager->GetNavigatorForTracking()->SetWorldVolume(UserDetectorConstruction->Construct());

  m_UserDetector          = UserDetectorConstruction;
  m_idmap                 = UserDetectorConstruction->idMap();
  
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

StatusCode G4GeometrySvc::queryInterface(const InterfaceID& riid, void** ppvInterface)  
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

idents::VolumeIdentifier G4GeometrySvc::getVolumeIdent(const G4VPhysicalVolume* volume)
{
  // Purpose and Method:  Given a volume return its identifier
  // Inputs:  Pointer to a G4VPhysicalVolume
  // Outputs:  The identifier associated with the volume
  // Dependencies: None
  // Restrictions None 
    return (*m_idmap)[volume];
}

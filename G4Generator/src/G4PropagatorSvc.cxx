#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "G4ParticlePropagator.h"

#include "GaudiKernel/Service.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/DataSvc.h"

#include "RunManager.h"
#include "DetectorConstruction.h"

#include <cassert>

/** 
* @class G4PropagatorSvc
*
* @brief A Gaudi Service for interfacing to the G4ParticlePropagator
*
* G4PropagatorSvc provides an interface to the world outside the G4Generator package
* which allows use of the G4ParticlePropagator for transporting tracks and their 
* error matrices throughout the GLAST detector.
* The current interface is taken from that used by RecoSvc, the Gismo equivalent
*
* @author Tracy Usher
*
*/

//class G4PropagatorSvc : virtual public DataSvc, virtual public IPropagatorSvc 
class G4PropagatorSvc : public IPropagatorSvc, public Service
{
public: 
    
    /// Return a pointer to the G4ParticlePropagator singleton object
    virtual IKalmanParticle* getPropagator() {return G4ParticlePropagator::instance();}

    virtual StatusCode initialize();

    virtual StatusCode finalize();
        
    /// queryInterface - for implementing a Service this is necessary
    virtual StatusCode queryInterface(const IID& riid, void** ppvUnknown);
    
protected:

    G4PropagatorSvc(const std::string& name, ISvcLocator* pSvcLocator);

    virtual ~G4PropagatorSvc();
  
private:
    // Allow SvcFactory to instantiate the service.
    friend class SvcFactory<G4PropagatorSvc>;

    //This is a pointer to the all important volume->idents map
    //obtained from the RunManager singleton
    //DetectorConstruction::IdMap* pIdMap;

    RunManager* m_RunManager;
};

static SvcFactory<G4PropagatorSvc> g4_factory;

const ISvcFactory& G4PropagatorSvcFactory = g4_factory;

//G4PropagatorSvc::G4PropagatorSvc(const std::string& name, ISvcLocator* pSvcLocator) :
//                 DataSvc(name, pSvcLocator)
G4PropagatorSvc::G4PropagatorSvc(const std::string& name, ISvcLocator* pSvcLocator) :
                 Service(name, pSvcLocator)
{
    return;
}

G4PropagatorSvc::~G4PropagatorSvc()
{
    //a bit tricky... may have been deleted already by G4Generator...
    if (RunManager::GetRunManager())
    {
        delete m_RunManager;
    }
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
    if (!(m_RunManager = RunManager::GetRunManager()))
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

        // The geant4 run manager
        m_RunManager = new RunManager(gsv,eventSvc);

        // Initialize Geant4
        m_RunManager->Initialize();

        log << MSG::DEBUG << "G4 RunManager ready" << endreq;
    }

    return sc;
}


StatusCode G4PropagatorSvc::finalize()
{
    // Purpose and Method:  Gaudi finalize routine. 
    // Inputs:  None
    // Outputs:  Gaudi StatusCode
    // Dependencies: None
    // Restrictions None 

    StatusCode sc = Service::finalize();

    return sc;
}


/// Query interface
StatusCode G4PropagatorSvc::queryInterface(const IID& riid, void** ppvInterface)  
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

// File and Version Information:
// $Header$

#define GlastApps_CreateEvent_CPP 

#include "CreateEvent.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"


#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/EventModel.h"
#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "GaudiKernel/ObjectVector.h"

#include "Event/MonteCarlo/McIntegratingHit.h"
#include "GlastSvc/Reco/IReconTool.h"

class MyRecon : public AlgTool, virtual public IReconTool {
 public:
     
   MyRecon( const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~MyRecon() { }
    
    /// implement to get something
    StatusCode get(const std::string& name, double&);
    
};


// Static factory for instantiation of algtool objects
static ToolFactory<MyRecon> afactory;

// Standard Constructor
MyRecon::MyRecon(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : AlgTool( type, name, parent ) {
    
    // Declare additional interface
    declareInterface<IReconTool>(this);
        
}


StatusCode MyRecon::get(const std::string& name, double& x)
{
    if(name=="test"){
        x = 99;
        return StatusCode::SUCCESS;
    }
    return StatusCode::FAILURE;
} 




static const AlgFactory<CreateEvent>  Factory;
const IAlgFactory& CreateEventFactory = Factory;

CreateEvent::CreateEvent(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator), m_detSvc(0) {
    
};

StatusCode CreateEvent::initialize() {
    // Purpose and Method: Test access to GlastDetSvc
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    // now try to find the GlastDevSvc service
    StatusCode sc = service("GlastDetSvc", m_detSvc, true);
    
    if (sc.isSuccess ()) {
        log << MSG::INFO << "Succeeded in accessing the GlastDetSvc!" << endreq;
        log << MSG::INFO << "testing constant access..." << endreq;
        double test;
        sc = m_detSvc->getNumericConstByName("junk", &test);
        if( sc.isFailure ()) log << MSG::INFO << "proper failure!" << endreq;
        sc = m_detSvc->getNumericConstByName("diodeX", &test);
        if( sc.isFailure ()) log << MSG::INFO << "diodeX not found!" << endreq;
        else log << MSG::INFO << "found constant diodeX = " << test << endreq;
        // Try out integer constant version
        int  myInt;
        sc = m_detSvc->getNumericConstByName("diodeX", &myInt);
        if( sc.isFailure ()) log << MSG::INFO << "proper failure!" << endreq;
        sc = m_detSvc->getNumericConstByName("xNum", &myInt);
        if (sc.isFailure()) log << MSG::INFO << "xNum not found!" << endreq;
        else log << MSG::INFO << "found constant xNum = " << myInt << endreq;
        
        // IDmap stuff
        idents::VolumeIdentifier prefix = m_detSvc->getIDPrefix();
        log << "Size of id prefix is " << prefix.size() << endreq;
        
        const HepTransform3D trans = m_detSvc->getTransform3DPrefix();
        const Hep3Vector vec = trans.getTranslation();
        log << "Prefix translation (x, y, z) is" << endreq;
        log << " (" << vec.x() << ", "
            << vec.y() << ", " << vec.z() << ")" << endreq;
        
    }else {
        log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
        return StatusCode::FAILURE;
    }
    
    // find the alg tool
    IToolSvc* tsvc  =0;
    sc = service( "ToolSvc", tsvc, true );
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to locate Tool Service" << endreq;
        return sc;
    }

    IReconTool* itool;
    sc = tsvc->retrieveTool("MyRecon", itool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find a test tool" << endreq;
        return sc;
    }

    double x=0;
    sc = itool->get("test",x);
    if( sc.isFailure() || x != 99) {
        log << MSG::ERROR << "Failed ReconTool test" << endreq;
    }
    return sc;
};

StatusCode CreateEvent::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << "execute" << endreq;
    
    //TODO: put something in here to get data???
    
    
    return sc;
    
};

StatusCode CreateEvent::finalize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize" << endreq;
    
    return StatusCode::SUCCESS;
}
;

StatusCode CreateEvent::testMcClass() {
    Event::McIntegratingHit* integratingHit = new Event::McIntegratingHit();
    Event::McParticle* mcParticle =  new Event::McParticle();
    
    
    return StatusCode::SUCCESS;
};




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
#include "CLHEP/Random/RandGauss.h"
#include "GaudiKernel/ObjectVector.h"

#include "Event/MonteCarlo/McIntegratingHit.h"
#include "GlastSvc/Reco/IReconTool.h"
#include "GlastSvc/GlastRandomSvc/RandomAccess.h"


/**
    @class MyRecon
    @brief example, test of mechanism
  */
class MyRecon : public AlgTool, virtual public IReconTool {
public:
    
    MyRecon( const std::string& type, 
            const std::string& name, 
            const IInterface* parent)
            :AlgTool(type,name,parent){
        
        // Declare additional interface
        declareInterface<IReconTool>(this);
        
    }
    virtual ~MyRecon() { }
    
    /// implement to get something
    StatusCode get(const std::string& name, double& x){
        if(name=="test"){
            x = 99;
            return StatusCode::SUCCESS;
        }
        return StatusCode::FAILURE;
    } 
    
};


// Static factory for instantiation of algtool objects
static ToolFactory<MyRecon> afactory;


//--------------------------------------------------------------------------
/**
    @class MyRandomTool
    @brief example, test of mechanism

    This makes a RandomAccess sub class that is instantiated in the local context,
    and then will allow the service to replace the pointer to the engine.
    
      TODO: put in some code to actually exercise setting the seed.
  */
class MyRandomTool : public RandomAccess {
public:
    MyRandomTool::MyRandomTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : RandomAccess( type, name, parent ) {
    
    }
};
static  ToolFactory<MyRandomTool>  randFactory;
//--------------------------------------------------------------------------
/**
    @class MyOtherRandomTool
    @brief example, more test of mechanism

    Add a second RandomAccess sub class (instantiated in the local context),
    to check that random engines actually get set at different addresses for
    different tools and different seeds get set for different tools
    
  */
class MyOtherRandomTool : public RandomAccess {
public:
    MyOtherRandomTool::MyOtherRandomTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : RandomAccess( type, name, parent ) {
    
    }
  
};
static  ToolFactory<MyOtherRandomTool>  otherRandFactory;


//--------------------------------------------------------------------------

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

        std::vector<idents::VolumeIdentifier> segs;
        m_detSvc->orderRibbonSegments(segs, 0, 1, false);
        log << "Top ribbon segments for ribbon 1, major axis in Y direction"
            << " (X-measuring) in order of increasing Y are:" << endreq;
        for (unsigned ix = 0; ix < segs.size(); ix++) {
          log <<  segs[ix].name() << endreq;
        }
        log << endreq;
        IGlastDetSvc::NamedId nid;
        nid = m_detSvc->getNamedId(segs[0]);
        log << "Named id corresponding to first segment is:" << endreq;
        for (unsigned ifld = 0; ifld < nid.size(); ifld++) {
          log << nid[ifld].first << "=" << nid[ifld].second << endreq;
        }
        log << endreq;
        
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
    log << MSG::DEBUG <<"Initial RandGauss status: " << RandGauss::getFlag() << endreq;
    // HMK Unused double test = RandGauss::shoot();
    log << MSG::DEBUG << "RandGauss status after a shoot: " << RandGauss::getFlag() << endreq;
    
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


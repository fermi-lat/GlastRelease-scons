// $Header$
#define GlastApps_CreateEvent_CPP 


//#define GlastApps_CreateEvent_CPP 

#include "CreateEvent.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/EventModel.h"
#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "GaudiKernel/ObjectVector.h"

#include "Event/MonteCarlo/McIntegratingHit.h"
                                
static const AlgFactory<CreateEvent>  Factory;
const IAlgFactory& CreateEventFactory = Factory;

//extern void GlastSvc_load();

/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.
CreateEvent::CreateEvent(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator), m_detSvc(0), m_irfLoadSvc(0) {
    //GlastSvc_load();
};


/*! The "functional" part of the class: For the EmptyAlgorithm example they do
nothing apart from print out info messages.
NB in the initialize method: you must explicitly initialize the base class
before using any services (message service, event data service etc.) otherwise 
the behaviour will be unpredictable at best.
*/
StatusCode CreateEvent::initialize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    // now try to find the GlastDevSvc service

    StatusCode sc = service("GlastDetSvc", m_detSvc);
    
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


    return StatusCode::SUCCESS;
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




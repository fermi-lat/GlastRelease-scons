/** @file BeamTransform.cxx
    @brief declartion, implementaion of the class BeamTransform

    $Header$
*/
// Gaudi system includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

// CLHEP for geometry
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"


// TDS class declarations: input data, and McParticle tree
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/MCEvent.h"


#include <cassert>
#include <sstream>

/** @class BeamTransform
@brief alg to transform the beam particle(s)
*/
class BeamTransform : public Algorithm {
public:
    BeamTransform(const std::string& name, ISvcLocator* pSvcLocator);
    /// set parameters and attach to various perhaps useful services.
    StatusCode initialize();
    /// process one event
    StatusCode execute();
    /// clean up
    StatusCode finalize();
    
private: 
    int m_count;
    DoubleProperty  m_transy, m_transz; // move in plane
    DoubleProperty  m_rotx, m_roty; //??

    void transform(Event::McParticle& mcp );
    HepRotation m_rotation;
    Hep3Vector m_translation;
};


static const AlgFactory<BeamTransform>  Factory;
const IAlgFactory& BeamTransformFactory = Factory;

BeamTransform::BeamTransform(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0)
{
    // declare properties with setProperties calls
    declareProperty("delta_y",  m_transy=0);
    declareProperty("delta_z",  m_transz=0);
    
}

StatusCode BeamTransform::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    // define the transformation matrix here.
    m_translation = Hep3Vector(0, m_transy, m_transz);
    return sc;
}
void BeamTransform::transform(Event::McParticle& mcp )
{
    
    // get the initial position and momentum
    HepLorentzVector initialMomentum = mcp.initialFourMomentum();
    HepPoint3D initialPosition = mcp.initialPosition();

    // transform them
    HepPoint3D newPosition = initialPosition + m_translation;

    HepLorentzVector newMomentum = initialMomentum;

    mcp.initialize(const_cast<Event::McParticle*>( &mcp.mother()), 
        mcp.particleProperty(), mcp.statusFlags(),
        newMomentum, newPosition,
        mcp.getProcess());
#if 0
    HepLorentzVector finalMomentum = mcp.finalFourMomentum();
    HepPoint3D finalPosition = mcp.finalPosition();
#endif

}

StatusCode BeamTransform::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    ++m_count;

    // Retrieving pointers from the TDS
    
    SmartDataPtr<Event::EventHeader>   header(eventSvc(),    EventModel::EventHeader);
    SmartDataPtr<Event::MCEvent>     mcheader(eventSvc(),    EventModel::MC::Event);
    SmartDataPtr<Event::McParticleCol> particles(eventSvc(), EventModel::MC::McParticleCol);

    double t = header->time();
    log << MSG::DEBUG << "Event time: " << t << endreq;;
  
    // loop over the monte carlo particle collection
    if (particles) {
        
        Event::McParticleCol::iterator piter;
        
        for (piter = particles->begin(); piter != particles->end(); piter++) {
            Event::McParticle& mcp = **piter;
            
            log << MSG::DEBUG ; log.stream() << mcp; log << endreq;
            
            // translate position
            transform(mcp);

            log << MSG::DEBUG << "after transform: "; log.stream() << mcp; log << endreq;
        }
    }

    return sc;
}

StatusCode BeamTransform::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "processed " << m_count << " events." << endreq;
    
    return sc;
}


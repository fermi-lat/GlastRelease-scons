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

// TU: Hack for CLHEP 1.9.2.2
typedef HepGeom::Point3D<double>  HepPoint3D;
typedef HepGeom::Vector3D<double> HepVector3D;

/** @class BeamTransform
@brief alg to transform the beam particle(s) into the detector frame

This involves a transformation of coordinates, and translation and rotation of the "x-y table."

Here's the latest word:

Beamline coordinate system:
     x is along the beam
     z is up
     y completes the right-handed system

CU coordinate system:
     z is opposite to the beam direction
     y is down
     x completes the right-handed system

The x-y table translations are specified in the beam system:
     m_transz is the up-down direction
     m_transy is the sideways direction

The pivot position is specified in the CU system
     m_pivot_height is in the z direction
     m_pivot_offset is in the x direction
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
    DoubleProperty m_pivot_height;    // z position of pivot in glast frame
    DoubleProperty m_pivot_offset;      // x position of pivot in glast frame
    DoubleProperty m_beam_plane; // z position of reporting plane in beam frame
    DoubleProperty m_beam_plane_glast;  // z position of reporting plane in glast frame
    DoubleProperty  m_transy, m_transz; // translation in x-y table plane
    DoubleProperty  m_angle; // rotation about table pivot (y-axis in instr.) (deg)
    double m_c, m_s;                    // cos, sin of rot

    void transform(Event::McParticle& mcp );
    CLHEP::Hep3Vector m_translation, m_pivot;
    CLHEP::HepRotationY m_rot;  // the table rotation
};

static const AlgFactory<BeamTransform>  Factory;
const IAlgFactory& BeamTransformFactory = Factory;

BeamTransform::BeamTransform(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0)
{
    // declare properties with setProperties calls
    // translations and rotation of x-y table
    declareProperty("vertical_translation",  m_transz=0);
    declareProperty("horizontal_translation",m_transy=0);
    declareProperty("table_rotation", m_angle=0);
    // z position of pivot in glast frame according to the Italians
    declareProperty("pivot_location", m_pivot_height=-130);
    // x position of pivot in glast frame according to the Italians
    declareProperty("pivot_offset", m_pivot_offset=45);
    // z position of reporting plane in beam frame
    declareProperty("beam_plane",     m_beam_plane=3300);
    // z coordinate of reporting plane in instrument frame
    declareProperty("beam_plane_glast" , m_beam_plane_glast=1000);
}

StatusCode BeamTransform::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    // define the transformation matrix here.
    m_translation = CLHEP::Hep3Vector(0, m_transy, m_transz);

    // this is the rotation matrix
    m_rot = CLHEP::HepRotationY(m_angle*M_PI/180);

    // location of pivot
    m_pivot = CLHEP::Hep3Vector( m_pivot_offset,0, m_pivot_height);

    return sc;
}
void BeamTransform::transform(Event::McParticle& mcp )
{
    // special treatment if it is a mother (code from Leon)
    if(&mcp.mother()!=&mcp) {
        const_cast<Event::McParticle*>(&mcp.mother())->removeDaughter(&mcp);
    }

    // get the initial position and momentum
    CLHEP::HepLorentzVector pbeam  = mcp.initialFourMomentum();
    CLHEP::Hep3Vector rbeam  = mcp.initialPosition();
    // ditto final
    CLHEP::HepLorentzVector pbeam1 = mcp.finalFourMomentum();
    CLHEP::Hep3Vector rbeam1 = mcp.finalPosition();

    // translate in beam frame
    rbeam  -= m_translation;
    rbeam1 -= m_translation;

    // convert to unrotated instrument coordinates
    CLHEP::Hep3Vector r (rbeam.y(),  -rbeam.z(), 
        -rbeam.x()  + m_beam_plane + m_beam_plane_glast);
    CLHEP::Hep3Vector r1(rbeam1.y(), -rbeam1.z(), 
        -rbeam1.x() + m_beam_plane + m_beam_plane_glast);

    CLHEP::HepLorentzVector p (pbeam.y(),  -pbeam.z(),  -pbeam.x(),  pbeam.e());
    CLHEP::HepLorentzVector p1(pbeam1.y(), -pbeam1.z(), -pbeam1.x(), pbeam1.e());

    mcp.initialize(const_cast<Event::McParticle*>( &mcp.mother()), 
        mcp.particleProperty(), mcp.statusFlags(),
        m_rot*p, 
        m_rot*(r-m_pivot) + m_pivot, 
        mcp.getProcess());

    mcp.finalize(
        m_rot*p1, 
        m_rot*(r1-m_pivot) + m_pivot);
}

StatusCode BeamTransform::execute(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    ++m_count;

    // Retrieving pointers from the TDS
    
    SmartDataPtr<Event::EventHeader>   header(eventSvc(),    EventModel::EventHeader);
    SmartDataPtr<Event::MCEvent>     mcheader(eventSvc(),    EventModel::MC::Event);
    SmartDataPtr<Event::McParticleCol> particles(eventSvc(), EventModel::MC::McParticleCol);

    double t = header->time();
    log << MSG::DEBUG;
    bool debug = log.isActive();
    //debug = true;
    if (debug) {log << "Event time: " << t;}
    log << endreq;;
  
    // loop over the monte carlo particle collection
    if (particles) {
        
        Event::McParticleCol::iterator piter;
        
        for (piter = particles->begin(); piter != particles->end(); piter++) {
            Event::McParticle& mcp = **piter;
            
            if (debug) {
                log << MSG::DEBUG;
                log.stream() << mcp;
                log << endreq;
            }
            
            // translate position
            transform(mcp);

            if (debug) {
                log << MSG::DEBUG << "after transform: "; 
                log.stream() << mcp; 
                log << endreq;
            }
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


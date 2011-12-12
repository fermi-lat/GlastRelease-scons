/** @file BeamTransform.cxx
    @brief declartion, implementaion of the class BeamTransform

    Transforms particle trajectories generated in the beamtest-2006 beam frame
    into the CU frame, taking into account the translations and rotations of the
    x-y table.

    @authors Toby Burnett, Leon Rochester

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

#include "geometry/Ray.h"

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
    /// z position of pivot in CU frame
    DoubleProperty m_pivot_height;    // z position of pivot in glast frame
    DoubleProperty m_pivot_offset;      // x position of pivot in glast frame
    DoubleProperty m_beam_plane; // z position of reporting plane in beam frame
    DoubleProperty m_beam_plane_glast;  // z position of reporting plane in glast frame
    DoubleProperty  m_transy, m_transz; // translation in x-y table plane
    DoubleProperty  m_angle; // rotation about table pivot (y-axis in instr.) (deg)
    DoubleProperty  m_tilt;  // rotation about x-axis in instr. (deg)
    DoubleArrayProperty m_point_on_beamline; // point in CU on central ray of beam
    double m_c, m_s;                    // cos, sin of rot angle
    double m_ct, m_st;                  // cos, sin of tilt

    void transform(Event::McParticle& mcp );
    CLHEP::Hep3Vector m_translation, m_pivot;
    CLHEP::HepRotationY m_rotY;  // the table rotation
    CLHEP::HepRotationX m_rotX;  // the table tilt
    CLHEP::HepRotation  m_rot;   // the full rotation
};

namespace {
    const double degToRad = M_PI/180.;
}

//static const AlgFactory<BeamTransform>  Factory;
//const IAlgFactory& BeamTransformFactory = Factory;
DECLARE_ALGORITHM_FACTORY(BeamTransform);

BeamTransform::BeamTransform(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0)
{
    // declare properties with setProperties calls
    // translations, rotation and tilt of x-y table
    //   distances in mm, angles in degrees
    // vertical translation of the table (along z(beam))
    declareProperty("vertical_translation",  m_transz=0);
    // horizontal translation of the table, (along y(beam))
    declareProperty("horizontal_translation",m_transy=0);
    // rotation of the table around the pivot
    declareProperty("table_rotation", m_angle=0);
    // tilt of the table around an axis along x(CU), going through the pivot
    declareProperty("table_tilt", m_tilt=0);
    // z position of pivot in glast frame
    declareProperty("pivot_location", m_pivot_height=-130);
    // x position of pivot in glast frame
    declareProperty("pivot_offset", m_pivot_offset=45);
    // z position of reporting plane in beam frame
    declareProperty("beam_plane",     m_beam_plane=3300);
    // z coordinate of reporting plane in instrument frame
    declareProperty("beam_plane_glast" , m_beam_plane_glast=1000);
    // point in CU along central ray of beam
    declareProperty("point_on_beamline", 
        m_point_on_beamline=std::vector<double>(3,-9999.));
}

StatusCode BeamTransform::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    const std::vector<double>& POBx = m_point_on_beamline.value( );
    Point POB = Point(POBx[0], POBx[1], POBx[2]);

        // this is the rotation matrix
    m_rotY = CLHEP::HepRotationY(m_angle*degToRad);
    m_rotX = CLHEP::HepRotationX(m_tilt*degToRad);
    m_rot  = m_rotX*m_rotY;

    // location of pivot
    m_pivot = CLHEP::Hep3Vector( m_pivot_offset,0, m_pivot_height);

    log << MSG::INFO << endreq << "*********************************************"
        << endreq << "*********************************************"
        << endreq << endreq;
    if(POB[0]!=-9999.) {
        log << MSG::INFO
            << "A point in the CU has been requested... " << endreq 
            << "The original x-y table coordinates will be over-written" << endreq;
        log << MSG::INFO << "CU point requested: ("
            << POB[0] << ", " << POB[1] << ", " << POB[2] << ")" 
            << ", angle: " << m_angle << " degrees" << endreq;
        log << MSG::INFO << "XY table Before: horizontal " << m_transy 
            << " mm, vertical " << m_transz << " mm" << endreq;
        m_transz = POB[1];
        if( m_angle==0&&m_tilt==0) {
            m_transy = -POB[0];
        } else {
            // I'm amazed that this code seems to work for angles of +-90.
            //    I guess that's the power of vector arithmetic!
            // That said, this seems a bit complicated...

            Point pivot(m_pivot.x(), m_pivot.y(), m_pivot.z());
            //make a ray orthogonal to the beam angle going thru the pivot

            double orthAngle = degToRad*(180. - m_angle);
            double tilt      = degToRad*m_tilt;
            Vector orthDir(cos(orthAngle), 0.0, sin(orthAngle));
            Ray orthRay(pivot, orthDir);

            // move by the offset of the pivot along this line
            // it's where the central ray enters for no translation
            Point point1 = orthRay.position(m_pivot.x());
           
            Vector dirBeam(0.0, 0.0, 1.0);
            dirBeam = m_rot*dirBeam;
            dirBeam = Vector(dirBeam.x(), 0.0, dirBeam.z());
            dirBeam.unit();

            Ray zeroRay(point1, dirBeam);

            // Here's the point you want the beam to go thru, projected onto the y(CU)=0 plane
            Point point2(POB[0], 0.0, POB[2]);
            Vector v1 = point2 - point1;
            // The distance along the beam central ray that gets you to the 
            //  doca for the beam point
            double lenB = v1*dirBeam;

            Point pointOnZero = zeroRay.position(lenB);
            Point projPOZ(pointOnZero.x(), 0.0, pointOnZero.z());
            Vector disp = point2-projPOZ;
            // This is how much you have to move
            double dist1 = disp.mag();
            // and this is the direction
            double sign = (disp*orthDir<0 ? -1.0 : 1.0);
            m_transy = sign*dist1;

            // correct vertical translation for the tilt
            m_transz = m_transz - (m_pivot.z()-POB[2])*tan(tilt);
        }
        log << MSG::INFO << "         After:  horizontal " << m_transy 
            << " mm, vertical " << m_transz << " mm" << endreq;
    } else {
        log << MSG::INFO << "XY table: horizontal " << m_transy 
            << " mm, vertical " << m_transz 
            << " mm, rotation: " << m_angle << " degrees" << endreq;
    }
    log << MSG::INFO << endreq << "*********************************************"
        << endreq << "*********************************************"
        << endreq << endreq;
 
    // define the transformation matrix here.
    m_translation = CLHEP::Hep3Vector(0, m_transy, m_transz);

    return sc;
}
void BeamTransform::transform(Event::McParticle& mcp )
{
    // get the positions and momenta
    // "m" means momentum, "c" means coordinates
    // "I" and "F" refer to initial and final quantities

    CLHEP::HepLorentzVector mBeamI = mcp.initialFourMomentum();
    CLHEP::Hep3Vector       cBeamI = mcp.initialPosition();
    // ditto final
    CLHEP::HepLorentzVector mBeamF = mcp.finalFourMomentum();
    CLHEP::Hep3Vector       cBeamF = mcp.finalPosition();

    // translate in beam frame
    cBeamI -= m_translation;
    cBeamF -= m_translation;

    // convert to unrotated instrument coordinates
    CLHEP::Hep3Vector cI(cBeamI.y(), -cBeamI.z(), 
        -cBeamI.x() + m_beam_plane + m_beam_plane_glast);
    CLHEP::Hep3Vector cF(cBeamF.y(), -cBeamF.z(), 
        -cBeamF.x() + m_beam_plane + m_beam_plane_glast);

    CLHEP::HepLorentzVector mI(mBeamI.y(), -mBeamI.z(), -mBeamI.x(), mBeamI.e());
    CLHEP::HepLorentzVector mF(mBeamF.y(), -mBeamF.z(), -mBeamF.x(), mBeamF.e());

    mcp.transform(m_rot*mI, m_rot*mF, 
        m_rot*(cI-m_pivot)+m_pivot, m_rot*(cF-m_pivot)+m_pivot);
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


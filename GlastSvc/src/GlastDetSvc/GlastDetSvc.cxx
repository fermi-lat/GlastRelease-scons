// File and Version Information:
// $Header$
// 
//  Original author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#include "GlastDetSvc.h"

#include <fstream>
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"

#include "GlastSvc/GlastDetSvc/IGeometry.h"
#include "GlastSvc/GlastDetSvc/IMedia.h"

#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Transform3D.h"

#include "SiliconPlaneGeometry.h"
#include "DMmanager.h"
#include "MediaVisitor.h"
#include "GeometryVisitor.h"

// define the class-id (unique) for the GlastDetSvc
// moved to Kernel! const IID&  IID_IGlastDetSvc  =  401;   // Unique to GLAST 

// define the class-id (unique) for the GlastIRFSvc
//const IID&  IID_IGlastIRFLoadSvc  =  402;   // Unique to GLAST


// declare the service factories for the GlastDetSvc
static SvcFactory<GlastDetSvc> a_factory;
const ISvcFactory& GlastDetSvcFactory = a_factory;

// external declarations of the service identifiers
// HSG - not necessary from the Gaudi example - extern const IID& IID_IGlastDetSvc;

GlastDetSvc::GlastDetSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc)
, m_instrument (0)
, m_xmlFile ("")
, m_iniFile ("")
, m_dm(0)
{
    
    // declare the properties
    declareProperty ("PersistencyFile", m_xmlFile);
    declareProperty ("IniFile", m_iniFile);

    // for DetModel visitor
    declareProperty("xmlfile",     m_xmlfile="-");
    declareProperty("topVolume",   m_topvol="LAT");
    declareProperty("visitorMode", m_visitorMode="propagate");
    
}

GlastDetSvc::~GlastDetSvc()  
{
    // this is  repetetive in case finalize was
    // not called.
    
    finalize();
}

StatusCode  GlastDetSvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_IGlastDetSvc == riid) {
        *ppvIF = dynamic_cast<IGlastDetSvc*> (this);
        return StatusCode::SUCCESS;
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
}


StatusCode GlastDetSvc::initialize () 
{
    
        // open the message log
    MsgStream log( msgSvc(), name() );

    StatusCode  status = StatusCode::SUCCESS;
    
    // call super-class
    Service::initialize ();
    
    // bind all of the properties for this service
    if ( (status = setProperties()).isFailure() ) {
        log << MSG::ERROR << "Failed to set properties" << endreq;
        return status;
    }
        
    // setup the detModel geometry, so can be visited below
    m_dm = new DMmanager;
    m_dm->init(m_xmlfile,m_visitorMode,m_topvol);
    log << MSG::INFO;
    m_dm->printSetup(log.stream());
    log << endreq;

    if( (status = SiliconPlaneGeometry::init(this)).isFailure() ) {
        log << MSG::WARNING << "Failed to initialize SiliconPlane Geometry " << endreq;
    }
    return StatusCode::SUCCESS;
}


StatusCode GlastDetSvc::finalize ()
{
    StatusCode  status = StatusCode::SUCCESS;
    
    
    // remove the detModel interface.
    delete m_dm;
    m_dm=0;
    return status;
}

const IID&  GlastDetSvc::type () const {
    return IID_IGlastDetSvc;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//            the new detModel interface calls
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GlastDetSvc::accept(IGeometry& geom)
{
    m_dm->accept(&GeometryVisitor(geom, m_visitorMode));
}

void GlastDetSvc::accept(IMedia& media)
{
  m_dm->accept(&MediaVisitor(media));
}

StatusCode GlastDetSvc::getNumericConstByName(std::string name, double* res)
{
  if (m_dm->getNumericConstByName(name, res))
    return StatusCode::SUCCESS;
  else return StatusCode::FAILURE;
}

StatusCode GlastDetSvc::getNumericConstByName(std::string name, int* res)
{
  if (m_dm->getNumericConstByName(name, res))
    return StatusCode::SUCCESS;
  else return StatusCode::FAILURE;
}

idents::VolumeIdentifier GlastDetSvc::getIDPrefix() 
{
  return m_dm->getIDPrefix();
}


StatusCode GlastDetSvc::getTransform3DByID(idents::VolumeIdentifier id,HepTransform3D* tr)
{
  if (m_dm->getTransform3DByID(id, tr))
    return StatusCode::SUCCESS;
  else return StatusCode::FAILURE;
}

const HepTransform3D& GlastDetSvc::getTransform3DPrefix()
{
  return m_dm->getTransform3DPrefix();
}

StatusCode  GlastDetSvc::getShapeByID(idents::VolumeIdentifier id,
                                      std::string* st, 
                                      std::vector<double>* dv)
{
  if (m_dm->getShapeByID(id, st, dv))
    return StatusCode::SUCCESS;
  else return StatusCode::FAILURE;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//            SiliconPlaneGeometry interface calls
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unsigned int GlastDetSvc::stripId (double x)
{
     return SiliconPlaneGeometry::stripId(x);
}

unsigned int GlastDetSvc::stripId (const HepPoint3D&  p)
{
     return SiliconPlaneGeometry::stripId(p.x());
}
/// (-) if non-active or gap (P is in local system)    
double GlastDetSvc::insideActiveArea (const HepPoint3D& p)
{
    return SiliconPlaneGeometry::insideActiveArea(p);
}

double GlastDetSvc::stripLocalX( double strip)
{
    return SiliconPlaneGeometry::localX(strip);
}

HepPoint3D GlastDetSvc::siPlaneCoord( const HepPoint3D &p, idents::VolumeIdentifier id)
{   
	return SiliconPlaneGeometry::siPlaneCoord( p, id);
}

HepPoint3D GlastDetSvc::getStripPosition(idents::VolumeIdentifier volId, double stripid)
{
    // Purpose: calculate the global position of a strip
    // Method:  apply the transform to the local position from SiliconPlaneGeometry
    // Inputs:  volId and strip
    // Output:  global position

    MsgStream log( msgSvc(), name() );

    HepTransform3D volTransform;
    StatusCode sc = getTransform3DByID(volId, &volTransform);
    if( sc.isFailure()) {
        log << MSG::INFO << "SiliconPlaneGeometry: Failed to obtain transform for id " 
           << volId.name() << endreq;
    }

    HepPoint3D localPos = SiliconPlaneGeometry::getLocalStripPosition(volId, stripid);

    return volTransform*localPos;

}

void GlastDetSvc::layerToTray (int layer, int view, int& tray, int& botTop)
{
    SiliconPlaneGeometry::layerToTray (layer, view, tray, botTop);
}

void GlastDetSvc::trayToLayer(int tray, int botTop, int& layer, int& view)
{
    SiliconPlaneGeometry::trayToLayer (layer, view, tray, botTop);
}

void GlastDetSvc::planeToLayer(int plane, int& layer, int& view)
{
    SiliconPlaneGeometry::planeToLayer (plane, layer, view);
}


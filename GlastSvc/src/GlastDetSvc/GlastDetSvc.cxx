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


// file scope for WARNING and ERROR
static MsgStream * s_log;

// external declarations of the service identifiers
// HSG - not necessary from the Gaudi example - extern const IID& IID_IGlastDetSvc;

// ------------------------------------------------
// Implementation of the GlastDetSvc class
// ------------------------------------------------
/// Standard Constructor
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
/// Standard Destructor
GlastDetSvc::~GlastDetSvc()  
{
    // this is  repetetive in case finalize was
    // not called.
    
    finalize();
}

// queryInterface
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


// initialize
StatusCode GlastDetSvc::initialize () 
{
    StatusCode  status = StatusCode::SUCCESS;
    
    // call super-class
    Service::initialize ();
    
    // bind all of the properties for this service
    status = setProperties ();
    
    
    // open the message log
    MsgStream log( msgSvc(), name() );
    s_log = & log;  // make available globally while executing the following
    
    
    // setup the detModel geometry, so can be visited below
    m_dm = new DMmanager;
    m_dm->init(m_xmlfile,m_visitorMode,m_topvol);
    log << MSG::INFO;
    m_dm->printSetup(log.stream());
    log << endreq;

	SiliconPlaneGeometry::init(this);

    return status;
}


// finalize
StatusCode GlastDetSvc::finalize ()
{
    StatusCode  status = StatusCode::SUCCESS;
    
    
    // remove the detModel interface.
    delete m_dm;
    m_dm=0;
    return status;
}

// access the type of this service
const IID&  GlastDetSvc::type () const {
    return IID_IGlastDetSvc;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//            the new detModel interface calls
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GlastDetSvc::accept(IGeometry& geom)
{
    m_dm->accept(&GeometryVisitor(geom));
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

StatusCode GlastDetSvc::getTransform3DByID(idents::VolumeIdentifier id,HepTransform3D* tr)
{
  if (m_dm->getTransform3DByID(id, tr))
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

    /// location of strip ix in local coords
double GlastDetSvc::stripLocalX ( unsigned int istrip)
{
    return SiliconPlaneGeometry::localX(istrip);
}

double GlastDetSvc::stripLocalXDouble( double strip)
{
    return SiliconPlaneGeometry::localXDouble(strip);
}

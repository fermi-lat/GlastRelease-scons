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
#include "xmlUtil/id/IdDict.h"

// define the class-id (unique) for the GlastDetSvc
// moved to Kernel! const InterfaceID&  IID_IGlastDetSvc  =  401;   // Unique to GLAST 

// define the class-id (unique) for the GlastIRFSvc
//const InterfaceID&  IID_IGlastIRFLoadSvc  =  402;   // Unique to GLAST


// declare the service factories for the GlastDetSvc
static SvcFactory<GlastDetSvc> a_factory;
const ISvcFactory& GlastDetSvcFactory = a_factory;

// external declarations of the service identifiers
// HSG - not necessary from the Gaudi example - extern const InterfaceID& IID_IGlastDetSvc;

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

StatusCode  GlastDetSvc::queryInterface (const InterfaceID& riid, void **ppvIF)
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

const InterfaceID&  GlastDetSvc::type () const {
    return IID_IGlastDetSvc;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//            the new detModel interface calls
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GlastDetSvc::accept(IGeometry& geom)
{
  GeometryVisitor gv(geom, m_visitorMode);
  m_dm->accept(&gv);
    //    m_dm->accept(&GeometryVisitor(geom, m_visitorMode));
}

void GlastDetSvc::accept(IMedia& media)
{
  MediaVisitor mv(media);
  m_dm->accept(&mv);
  //  m_dm->accept(&MediaVisitor(media));
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

IGlastDetSvc::NamedId 
GlastDetSvc::getNamedId(const idents::VolumeIdentifier& vid) {
  // first create an xmlUtil::Identifier - vector of unsigned - from
  // the volume identifier.
  xmlUtil::Identifier identifier;
  int iField;
  for (iField = 0; iField < vid.size(); iField++) {
    identifier.append(vid[iField]);
  }

  // then..

  xmlUtil::IdDict* pDict = m_dm->getIdDictionary();
  xmlUtil::NamedId *pId = pDict->getNamedId(identifier);
  IGlastDetSvc::NamedId nid = pId->fieldVector();
  delete pId;
  return nid;
}

void 
GlastDetSvc::orderRibbonSegments(std::vector<idents::VolumeIdentifier>& segs,
                                 unsigned face, unsigned ribbonNumber, 
                                 bool xOrient, bool increasing) 
{
  m_dm->orderRibbonSegments(segs, face, ribbonNumber, xOrient, increasing);
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

double GlastDetSvc::insideActiveLocalX( const HepPoint3D& p)
{
    return SiliconPlaneGeometry::insideActiveArea1D(p.x(), SiliconPlaneGeometry::ladder_gap());
}

double GlastDetSvc::insideActiveLocalY( const HepPoint3D& p)
{
    return SiliconPlaneGeometry::insideActiveArea1D(p.y(), SiliconPlaneGeometry::ssd_gap());
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

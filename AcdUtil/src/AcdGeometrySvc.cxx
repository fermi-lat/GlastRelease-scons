
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "src/AcdGeometrySvc.h"

#include "idents/TowerId.h"

#include <iostream>
#include <algorithm>

static const SvcFactory<AcdGeometrySvc> s_factory;
const ISvcFactory& AcdGeometrySvcFactory = s_factory;

//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

AcdGeometrySvc::AcdGeometrySvc(const std::string& name, 
                               ISvcLocator* pSvcLocator) : Service(name, pSvcLocator)
{   
    clear();
    return; 
}

StatusCode AcdGeometrySvc::initialize()
{
    // Purpose: load up constants from GlastDetSvc and do some calcs
    // Inputs:  none
    // Output:  AcdGeometrySvc statics initialized

    StatusCode sc = StatusCode::SUCCESS;

    Service::initialize();
    setProperties();
    MsgStream log(msgSvc(), name());

    m_glastDetSvc = 0;
    if (service("GlastDetSvc", m_glastDetSvc, true).isFailure()) {
        log << MSG::ERROR << "GlastDetSvc not found" << endreq;
        return StatusCode::FAILURE;
    }

    sc = getConstants();

    sc = getDetectorListFromGeometry();

    log << MSG::INFO << "AcdGeometrySvc successfully initialized" << endreq;
    return StatusCode::SUCCESS;

}

void AcdGeometrySvc::clear() {
    m_numXtowers = 0; m_numYtowers = 0;
    m_numTiles = 0;  m_numRibbons = 0;
}

StatusCode AcdGeometrySvc::finalize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "AcdGeometrySvc finalize called" << endreq;
    return StatusCode::SUCCESS;
}


StatusCode  AcdGeometrySvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_IAcdGeometrySvc == riid) {
        *ppvIF = dynamic_cast<IAcdGeometrySvc*> (this);
        return StatusCode::SUCCESS;
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
}

// access the type of this service

const IID&  AcdGeometrySvc::type () const {
    return IID_IAcdGeometrySvc;
}


StatusCode AcdGeometrySvc::getConstants()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    if (m_glastDetSvc->getNumericConstByName("xNum", &m_numXtowers).isSuccess() &&
        m_glastDetSvc->getNumericConstByName("yNum", &m_numYtowers).isSuccess()        ) 
     {
       sc = StatusCode::SUCCESS;
     } else {
        log << MSG::ERROR << "Failed to get geometry constants" << endreq;
        return StatusCode::FAILURE;
     }
    return sc;
}

StatusCode AcdGeometrySvc::getDetectorListFromGeometry() {

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // get the list of layers, to be used to add noise to otherwise empty layers
    m_detectorCol.setPrefix(m_glastDetSvc->getIDPrefix());

    // Find all the ACD detectors in our geometry
    m_glastDetSvc->accept(m_detectorCol);
    if (m_detectorCol.size() > 0)
        log << MSG::INFO << "Located  "<< m_detectorCol.size() 
            << " ACD volumes, ids from " << m_detectorCol.front().name() 
            << " to " << m_detectorCol.back().name() << endreq;
    else
        log << MSG::INFO << "No ACD detectors in this geometry" << endreq;

   for(AcdUtil::AcdDetectorList::const_iterator it=m_detectorCol.begin(); 
                                                it!=m_detectorCol.end(); ++it)
   {
       idents::VolumeIdentifier volId = *it;
       idents::AcdId detectorId(volId);
       if (detectorId.ribbon()) ++m_numRibbons;
       else if (detectorId.tile()) ++m_numTiles;
   }
   log << MSG::INFO << "Number of Tiles:  " << m_numTiles
       << " Number of Ribbons: " << m_numRibbons << endreq;

    return sc;
}


bool AcdGeometrySvc::findDetector(const idents::VolumeIdentifier &volId) const {
    std::vector<idents::VolumeIdentifier>::const_iterator foundIt;
    foundIt = find(m_detectorCol.begin(), m_detectorCol.end(), volId);
    return (foundIt != m_detectorCol.end());
}


StatusCode AcdGeometrySvc::getDimensions(
                           const idents::VolumeIdentifier &volId,
                           std::vector<double> &dims, HepPoint3D &xT) const {
    MsgStream   log( msgSvc(), name() );
    std::string str;
    StatusCode sc = StatusCode::SUCCESS;
    sc = m_glastDetSvc->getShapeByID(volId, &str, &dims);
    if ( sc.isFailure() ) {
        log << MSG::WARNING << "Failed to retrieve Shape by Id" << endreq;
        return sc;
    }
    HepTransform3D transform;
    sc = m_glastDetSvc->getTransform3DByID(volId, &transform);
    if (sc.isFailure() ) {
        log << MSG::WARNING << "Failed to get transformation" << endreq;
        return sc;
    }
    HepPoint3D center(0., 0., 0.);
    xT = transform * center;
    return sc;

}


StatusCode AcdGeometrySvc::getDetectorDimensions(
                           const idents::VolumeIdentifier &volId,
                           std::vector<double> &dims, HepPoint3D &xT) const {
    // Purpose and Method: Retrieve the dimensions - and take note of rotation
    //  of the side tiles, and return X, Y, Z as we would expect

    StatusCode sc;
    sc = getDimensions(volId, dims, xT);
    if (sc.isFailure()) return sc;
    idents::AcdId id(volId);
    double dX = dims[0];
    double dY = dims[1];
    if (id.face() == 1 || id.face() == 3) {
        dims[0] = dY;
        dims[1] = dX;
    }
    return sc;
}


StatusCode AcdGeometrySvc::getCorners(const std::vector<double> &dim,
                                const HepPoint3D &center, HepPoint3D *corner) {
    StatusCode sc = StatusCode::SUCCESS;

    unsigned int iCorner;
    // Ignore short dimension - only interested in 4 corners

    if ((dim[0] < dim[1]) && (dim[0] < dim[2])) { // X smallest -  X Face
        for(iCorner = 0; iCorner<4; iCorner++) {
            corner[iCorner].setX(center.x());
            corner[iCorner].setY(center.y() +
                     ( (iCorner < 2) ? -1 : 1) * dim[1]*0.5);
            corner[iCorner].setZ(center.z() +
                     ( ((iCorner == 0) || (iCorner==3)) ? -1 : 1) * dim[2]*0.5);
        }

    } else if ((dim[1] < dim[0]) && (dim[1] < dim[2])) {
        for(iCorner = 0; iCorner<4; iCorner++) {
            corner[iCorner].setY(center.y());
            corner[iCorner].setX(center.x() +
                     ( (iCorner < 2) ? -1 : 1) * dim[0]*0.5);
            corner[iCorner].setZ(center.z() +
                     ( ((iCorner == 0) || (iCorner==3)) ? -1 : 1) * dim[2]*0.5);
        }

    } else {
        for(iCorner = 0; iCorner<4; iCorner++) {
            corner[iCorner].setZ(center.z());
            corner[iCorner].setX(center.x() +
                     ( (iCorner < 2) ? -1 : 1) * dim[0]*0.5);
            corner[iCorner].setY(center.y() +
                     ( ((iCorner == 0) || (iCorner==3)) ? -1 : 1) * dim[1]*0.5);
        }

    }

    return sc;
}



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
    unsigned int i;
    for (i=0; i<4; i++) {
        m_cornerGapStartPoint[i] = Point(0.0, 0.0, 0.0);
        m_cornerGapEndPoint[i] = Point(0.0, 0.0, 0.0);
        m_cornerGapVec[i] = Vector(0.0, 0.0, 0.0);
    }
}

StatusCode AcdGeometrySvc::finalize()
{
//    MsgStream log(msgSvc(), name());
//    log << MSG::INFO << "AcdGeometrySvc finalize called" << endreq;
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

    if (m_glastDetSvc->getNumericConstByName("xNum", &m_numXtowers).isSuccess() &&
        m_glastDetSvc->getNumericConstByName("yNum", &m_numYtowers).isSuccess()        ) 
     {
       sc = StatusCode::SUCCESS;
     } else {
        MsgStream log(msgSvc(), name());
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
    std::string str;
    StatusCode sc = StatusCode::SUCCESS;
    sc = m_glastDetSvc->getShapeByID(volId, &str, &dims);
    if ( sc.isFailure() ) {
        MsgStream   log( msgSvc(), name() );
        log << MSG::WARNING << "Failed to retrieve Shape by Id" << endreq;
        return sc;
    }
    HepTransform3D transform;
    sc = m_glastDetSvc->getTransform3DByID(volId, &transform);
    if (sc.isFailure() ) {
        MsgStream   log( msgSvc(), name() );
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

const Ray AcdGeometrySvc::getCornerGapRay(unsigned int index) const {
    if (index > 4) {
        MsgStream  log( msgSvc(), name() );
        log << MSG::WARNING << "Requested Ray index out of bounds" << endreq;
        return Ray(Point(0,0,0), Vector(0,0,0));
    }
    
   Ray r(m_cornerGapStartPoint[index], m_cornerGapVec[index]);
   r.setArcLength(m_cornerGapVec[index].mag());
   return r;
}

StatusCode AcdGeometrySvc::findCornerGaps( ) {
    MsgStream   log( msgSvc(), name() );
    StatusCode sc = StatusCode::SUCCESS;

    // Form Corner Gap Vectors
    // Shameless hack..we know AcdIds for the corners:
    // 100, 200, 204, 300, 304, 404, 400, 104


    // Form the AcdIds of interest
    idents::AcdId id100(0, 1, 0, 0);
    idents::AcdId id104(0, 1, 0, 4);
    idents::AcdId id200(0, 2, 0, 0);
    idents::AcdId id204(0, 2, 0, 4);
    idents::AcdId id300(0, 3, 0, 0);
    idents::AcdId id304(0, 3, 0, 4);
    idents::AcdId id400(0, 4, 0, 0);
    idents::AcdId id404(0, 4, 0, 4);

    // Get one top tile and one bottom to figure out z top and bottom
    idents::AcdId id040(0, 0, 4, 0);
    idents::AcdId id430(0, 4, 3, 0);

    // Now get the VolumeIds
    idents::VolumeIdentifier volId100 = id100.volId();
    idents::VolumeIdentifier volId104 = id104.volId();
    idents::VolumeIdentifier volId200 = id200.volId();
    idents::VolumeIdentifier volId204 = id204.volId();
    idents::VolumeIdentifier volId300 = id300.volId();
    idents::VolumeIdentifier volId304 = id304.volId();
    idents::VolumeIdentifier volId400 = id400.volId();
    idents::VolumeIdentifier volId404 = id404.volId();
    idents::VolumeIdentifier volId040 = id040.volId();
    idents::VolumeIdentifier volId430 = id430.volId();

   // Check to see if these volIds exist in the geometry
    if (!findDetector(volId100)) sc = StatusCode::FAILURE;
    if (!findDetector(volId104)) sc = StatusCode::FAILURE;
    if (!findDetector(volId200)) sc = StatusCode::FAILURE;
    if (!findDetector(volId204)) sc = StatusCode::FAILURE;
    if (!findDetector(volId300)) sc = StatusCode::FAILURE;
    if (!findDetector(volId304)) sc = StatusCode::FAILURE;
    if (!findDetector(volId400)) sc = StatusCode::FAILURE;
    if (!findDetector(volId404)) sc = StatusCode::FAILURE;
    if (!findDetector(volId040)) sc = StatusCode::FAILURE;
    if (!findDetector(volId430)) sc = StatusCode::FAILURE;

    // If we fail to find any of the tiles, we skip calculation of the corner
    // gap variable
    if (sc.isFailure()) {
        log << MSG::WARNING << "Non-flight ACD Geometry will not calculate"
            << " corner gap rays" << endreq;
        return sc;
    }

    // Determine the extent of all corner gap rays in Z, using the a top corner
    // tile and a bottom side row tile.
    std::vector<double> dim040, dim430;
    HepPoint3D centerZ, cornerZ[4];
    // Using new function that automatically handles flipping X and Y for
    // Faces 1 and 3
    sc = getDetectorDimensions(volId040, dim040, centerZ);
    getCorners(dim040, centerZ, cornerZ);
    double zMax = cornerZ[0].z();
    sc = getDetectorDimensions(volId430, dim430, centerZ);
    getCorners(dim430, centerZ, cornerZ);
    double zMin = cornerZ[0].z();

    // Construct Ray for 1st Corner Gap
    std::vector<double> dim100, dim200;
    HepPoint3D center100, center200;
    HepPoint3D corner100[4], corner200[4];
    sc = getDetectorDimensions(volId100, dim100, center100);
    getCorners(dim100, center100, corner100);
    // use corner with minX and maxZ which should be index 1

    sc = getDetectorDimensions(volId200, dim200, center200);
    getCorners(dim200, center200, corner200);
    // use corner with minX and maxZ which should be index 1
    // Construct the Ray with the x and y position in the middle of the gap
    double xDiff = fabs(corner200[1].x() - corner100[1].x())*0.5;
    double yDiff = fabs(corner200[1].y() - corner100[1].y())*0.5;

    // Create the start and end of the ray using the position half-way between
    // the two corners of the adjacent tiles
    double xMin = corner100[1].x();
    double yMin = corner100[1].y();
    if (corner200[1].x() < corner100[1].x()) xMin = corner200[1].x();
    if (corner200[1].y() < corner100[1].y()) yMin = corner200[1].y();
    m_cornerGapStartPoint[0] = Point(xMin + xDiff, yMin + yDiff, zMin);
    m_cornerGapEndPoint[0] = Point(xMin + xDiff, yMin + yDiff, zMax);
    m_cornerGapVec[0] = Vector(m_cornerGapEndPoint[0] - m_cornerGapStartPoint[0]);
    //m_cornerGapRay[0] = Ray(m_cornerGapStartPoint[0], cornerGapVec[0]);
    log << MSG::DEBUG << "Corner1Ray:  Pos:(" << m_cornerGapStartPoint[0].x() 
        << ", " << m_cornerGapStartPoint[0].y() << ", " 
        << m_cornerGapStartPoint[0].z() << ") Vec:( "
        << m_cornerGapVec[0].x() << ", " << m_cornerGapVec[0].y() << ", " 
        << m_cornerGapVec[0].z() << ") mag: " << m_cornerGapVec[0].mag() << endreq;

    // Construct Ray for 2nd Corner Gap
    std::vector<double> dim300, dim204;
    HepPoint3D center300, center204;
    HepPoint3D corner300[4], corner204[4];
    // Using new function that automatically handles flipping X and Y for
    // Faces 1 and 3
    sc = getDetectorDimensions(volId300, dim300, center300);
    getCorners(dim300, center300, corner300);
    // use corner with minY and maxZ which should be index 1

    sc = getDetectorDimensions(volId204, dim204, center204);
    getCorners(dim204, center204, corner204);
    // use corner with maxX and maxZ which should be index 2
    // Construct the Ray with the x and y position in the middle of the gap
    xDiff = fabs(corner204[2].x() - corner300[1].x())*0.5;
    yDiff = fabs(corner204[2].y() - corner300[1].y())*0.5;

    // Create the start and end of the ray using the position half-way between
    // the two corners of the adjacent tiles
    xMin = corner300[1].x();
    yMin = corner300[1].y();
    if (corner204[2].x() < corner300[1].x()) xMin = corner204[2].x();
    if (corner204[2].y() < corner300[1].y()) yMin = corner204[2].y();
    m_cornerGapStartPoint[1] = Point(xMin + xDiff, yMin + yDiff, zMin);
    m_cornerGapEndPoint[1] = Point(xMin + xDiff, yMin + yDiff, zMax);
    m_cornerGapVec[1] = Vector(m_cornerGapEndPoint[1] - m_cornerGapStartPoint[1]);
    log << MSG::DEBUG << "Corner2Ray:  Pos:(" << m_cornerGapStartPoint[1].x() 
        << ", " << m_cornerGapStartPoint[1].y() << ", " 
        << m_cornerGapStartPoint[1].z() << ") Vec:( "
        << m_cornerGapVec[1].x() << ", " << m_cornerGapVec[1].y() << ", " 
        << m_cornerGapVec[1].z() << ") mag: " 
        << m_cornerGapVec[1].mag() << endreq;


   // Construct Ray for 3rd Corner Gap
    std::vector<double> dim404, dim304;
    HepPoint3D center404, center304;
    HepPoint3D corner404[4], corner304[4];
    sc = getDetectorDimensions(volId304, dim304, center304);
    getCorners(dim304, center304, corner304);
    // use corner with maxY and maxZ which should be index 2

    sc = getDetectorDimensions(volId404, dim404, center404);
    getCorners(dim404, center404, corner404);
    // use corner with maxX and maxZ which should be index 2
    // Construct the Ray with the x and y position in the middle of the gap
    xDiff = fabs(corner304[2].x() - corner404[2].x())*0.5;
    yDiff = fabs(corner304[2].y() - corner404[2].y())*0.5;

    // Create the start and end of the ray using the position half-way between
    // the two corners of the adjacent tiles
    xMin = corner404[2].x();
    yMin = corner404[2].y();
    if (corner304[2].x() < corner404[2].x()) xMin = corner304[2].x();
    if (corner304[2].y() < corner404[2].y()) yMin = corner304[2].y();
    m_cornerGapStartPoint[2] = Point(xMin + xDiff, yMin + yDiff, zMin);
    m_cornerGapEndPoint[2] = Point(xMin + xDiff, yMin + yDiff, zMax);
    m_cornerGapVec[2] = Vector(m_cornerGapEndPoint[2] - m_cornerGapStartPoint[2]);
    log << MSG::DEBUG << "Corner3Ray:  Pos:(" << m_cornerGapStartPoint[2].x() 
        << ", " << m_cornerGapStartPoint[2].y() << ", " 
        << m_cornerGapStartPoint[2].z() << ") Vec:( "
        << m_cornerGapVec[2].x() << ", " << m_cornerGapVec[2].y() << ", " 
        << m_cornerGapVec[2].z() << ") mag: " 
        << m_cornerGapVec[2].mag() << endreq;


   // Construct Ray for 4th Corner Gap
    std::vector<double> dim104, dim400;
    HepPoint3D center104, center400;
    HepPoint3D corner104[4], corner400[4];

    // Using new AcdGeometrySvc version to get dimensions
    sc = getDetectorDimensions(volId104, dim104, center104);
    // The corners are returned in order (-,-), (-,+), (+,+), (+,-)
    // where third dimension is the one we ignore, since it is associated with
    // tile thickness.
    getCorners(dim104, center104, corner104);
    // use corner with max Y and max Z, which should be index 2

    sc = getDetectorDimensions(volId400, dim400, center400);
    getCorners(dim400, center400, corner400);
    // Use corner with min X and max Z which should be index 1

    // Construct the Ray with the x and y position in the middle of the gap
    xDiff = fabs(corner104[2].x() - corner400[1].x())*0.5;
    yDiff = fabs(corner104[2].y() - corner400[1].y())*0.5;

    // Create the start and end of the ray using the position half-way between
    // the two corners of the adjacent tiles
    xMin = corner104[2].x();
    yMin = corner104[2].y();
    if (corner400[1].x() < corner104[2].x()) xMin = corner400[1].x();
    if (corner400[1].y() < corner104[2].y()) yMin = corner400[1].y();
    m_cornerGapStartPoint[3] = Point(xMin + xDiff, yMin + yDiff, zMin);
    m_cornerGapEndPoint[3] = Point(xMin + xDiff, yMin + yDiff, zMax);

    m_cornerGapVec[3] = Vector(m_cornerGapEndPoint[3] - m_cornerGapStartPoint[3]);
    log << MSG::DEBUG << "Corner4Ray:  Pos:(" << m_cornerGapStartPoint[3].x() 
        << ", " << m_cornerGapStartPoint[3].y() << ", " 
        << m_cornerGapStartPoint[3].z() << ") Vec:( "
        << m_cornerGapVec[3].x() << ", " << m_cornerGapVec[3].y() << ", " 
        << m_cornerGapVec[3].z() << ") mag: " 
        << m_cornerGapVec[3].mag() << endreq;

    return sc;
}



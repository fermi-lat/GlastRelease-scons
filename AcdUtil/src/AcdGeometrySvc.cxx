
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "src/AcdGeometrySvc.h"

#include "CLHEP/Geometry/Transform3D.h"

#include "idents/TowerId.h"

#include <iostream>
#include <algorithm>

static const SvcFactory<AcdGeometrySvc> s_factory;
const ISvcFactory& AcdGeometrySvcFactory = s_factory;

//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

AcdGeometrySvc::AcdGeometrySvc(const std::string& name, 
                               ISvcLocator* pSvcLocator) 
                              : Service(name, pSvcLocator)
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
    m_acdId_volCount.clear();
}

StatusCode AcdGeometrySvc::finalize()
{
//    MsgStream log(msgSvc(), name());
//    log << MSG::INFO << "AcdGeometrySvc finalize called" << endreq;
    return StatusCode::SUCCESS;
}


StatusCode AcdGeometrySvc::queryInterface(const InterfaceID& riid, void **ppvIF)
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
const InterfaceID&  AcdGeometrySvc::type () const {
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

       // Keep a count of volumes associated with each AcdId
       if (m_acdId_volCount.find(detectorId) != m_acdId_volCount.end()) {
           m_acdId_volCount[detectorId] += 1;
       } else {
           m_acdId_volCount[detectorId] = 1;
       }
   }


   log << MSG::DEBUG << "Found " << m_acdId_volCount.size() << " volumes"
       << endreq;
   std::map<idents::AcdId, int>::const_iterator volColIt;
   for (volColIt = m_acdId_volCount.begin(); volColIt != m_acdId_volCount.end();
        volColIt++) {
       idents::AcdId id = volColIt->first;
       if (id.ribbon()) ++m_numRibbons;
       else if (id.tile()) ++m_numTiles;
       log << MSG::DEBUG << id.id() << " has " << volColIt->second 
           << " volumes associated with it" << endreq;
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
        log << MSG::WARNING << "Failed to retrieve Shape by Id: " 
            << volId.name() << endreq;
        return sc;
    }
    HepGeom::Transform3D transform;
    sc = m_glastDetSvc->getTransform3DByID(volId, &transform);
    if (sc.isFailure() ) {
        MsgStream   log( msgSvc(), name() );
        log << MSG::WARNING << "Failed to get transformation: " 
            << volId.name() << endreq;
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

bool AcdGeometrySvc::fillRibbonRays(const idents::AcdId& id, std::vector<Ray>& minusSideRays, std::vector<Ray>& topRays,
                                    std::vector<Ray>& plusSideRays, bool increasing)
{
    // Purpose and Method:  Fill the three supplied vector of Rays.  The Rays are constructed from the ribbon segments
    //    associated with AcdId id.  

    minusSideRays.clear();
    topRays.clear();
    plusSideRays.clear();

    typedef enum {
        ribbonX = 5,
        ribbonY = 6
    } ribbonOrient;

    unsigned int faceArr[6] = {0, 1, 3, 0, 2, 4}; // x Faces 0, 1 (-X), 3 (+X), y Faces 0, 2 (-y), 4(+y)
    // do top faces first too add segments 6 and 7 to the right end of the minus and plus vectors
    // when increasing, it means we can just push_back the rays as normal, and when decreasing we can also use push_back

    bool retVal = true;
    if (!id.ribbon()) return false;

    unsigned int ribbonNum = id.ribbonNum();
    unsigned int ribbonOrientation = id.ribbonOrientation();

    // Use orientation of ribbon to determine xOrient and the starting index into our array of face ids
    bool xOrient;
    unsigned int startInd = 0;
    if (ribbonOrientation == ribbonX) 
        xOrient = true;
    else {
        startInd = 3;
        xOrient = false;
    }

    // Loops over set of faces
    unsigned int iFace;
    for (iFace = startInd; iFace < (startInd+3); iFace++) {

        // Get the ribbon segments from the glastDetSvc
        std::vector<idents::VolumeIdentifier> ribbonSegmentVolIds;
        m_glastDetSvc->orderRibbonSegments(ribbonSegmentVolIds,
            faceArr[iFace], ribbonNum, xOrient, increasing);

        // Determine the dimension to use to contruct the rays
        std::vector<idents::VolumeIdentifier>::const_iterator it;

        if ( ribbonSegmentVolIds.size() == 0 ) {
            MsgStream   log( msgSvc(), name() );
            log << MSG::INFO << "No ribbon segments found for AcdId: " << id.id() << " with xOrientation "
                << xOrient << " Face: " << faceArr[iFace] << endreq;
        }
    
        std::vector<double> dims;
        // Iterate over ribbon segments and construct rays
        for (it = ribbonSegmentVolIds.begin(); it != ribbonSegmentVolIds.end(); it++) {
            dims.clear();
            HepPoint3D center;
            unsigned int dimInd;

            if ( faceArr[iFace] == 0 ) { // On the top, interested in segments 1,2,3,4,5
                if ( (*it)[5] <= 5 ) {
                    dimInd = (xOrient) ? 0 : 1;
                } else
                    continue;
            } else {
                dimInd = 2;
                if ( ( (*it)[5] > 4) && ( (*it)[5] != 9) ) // On the sides we are interested in segments 1,2,3,4,9
                    continue;
            }

            
            getDimensions(*it, dims, center); 
            double cen[3] = {center.x(), center.y(), center.z()};
            Point startPos(cen[0], cen[1], cen[2]);
            Point endPos(cen[0], cen[1], cen[2]);
            if (increasing) {
               startPos(dimInd) = cen[dimInd] - dims[dimInd]/2.;
               endPos(dimInd) = cen[dimInd] + dims[dimInd]/2.;
            } else {
                startPos(dimInd) = cen[dimInd] + dims[dimInd]/2.;
                endPos(dimInd) = cen[dimInd] - dims[dimInd]/2.;
            }

            // Construct Rays, forming vector from end and starting position along ribbon segment
            Vector rayVec = endPos - startPos;
            Ray r(startPos, rayVec);
            r.setArcLength(rayVec.mag());
            MsgStream   log( msgSvc(), name() );
            log << MSG::DEBUG << "dimInd: " << dimInd << endreq;
            log << MSG::DEBUG << "cen: ( " << cen[0] << ", " << cen[1] << ", " << cen[2] << ")" << endreq;
            log << MSG::DEBUG << "dims: " << dims[0] << ", " << dims[1] << ", " << dims[2] << endreq;
            log << MSG::DEBUG << "startPos: (" << startPos.x() << ", " << startPos.y() << ", " << startPos.z() << ")" 
                << " endPos: ( " << endPos.x() << ", " << endPos.y() << ", " << endPos.z() << ")" << endreq;
 
            if (faceArr[iFace] == 0) {             // Face 0, segment# <= 5
                topRays.push_back(r);
            } else if ((faceArr[iFace] == 1) || (faceArr[iFace] == 2))  // Faces 1, 2
                minusSideRays.push_back(r);
            else                                                        // Faces 3,4
                plusSideRays.push_back(r);
    
        } // end ribbonSegment for

    } // end face for

    return retVal;
}


/// Given an AcdId for a ribbon, provide the transformation to the center of each set of ribbon segments
bool AcdGeometrySvc::fillRibbonTransforms(const idents::AcdId& id,
					  HepTransform3D& minusSideTransform,
					  HepTransform3D& topTransform,
					  HepTransform3D& plusSideTransform){

    typedef enum {
        ribbonX = 5,
        ribbonY = 6
    } ribbonOrient;

    unsigned int faceArr[6] = {0, 1, 3, 0, 2, 4}; // x Faces 0, 1 (-X), 3 (+X), y Faces 0, 2 (-y), 4(+y)
    // do top faces first too add segments 6 and 7 to the right end of the minus and plus vectors
    // when increasing, it means we can just push_back the rays as normal, and when decreasing we can also use push_back

    bool retVal = true;
    if (!id.ribbon()) return false;

    unsigned int ribbonNum = id.ribbonNum();
    unsigned int ribbonOrientation = id.ribbonOrientation();

    // Use orientation of ribbon to determine xOrient and the starting index into our array of face ids
    bool xOrient;
    unsigned int startInd = 0;
    if (ribbonOrientation == ribbonX) 
        xOrient = true;
    else {
        startInd = 3;
        xOrient = false;
    }

    // Loops over set of faces
    unsigned int iFace;
    for (iFace = startInd; iFace < (startInd+3); iFace++) {

        // Get the ribbon segments from the glastDetSvc
        std::vector<idents::VolumeIdentifier> ribbonSegmentVolIds;
        m_glastDetSvc->orderRibbonSegments(ribbonSegmentVolIds,
					   faceArr[iFace], ribbonNum, xOrient, true);

        if ( ribbonSegmentVolIds.size() == 0 ) {
            MsgStream   log( msgSvc(), name() );
            log << MSG::INFO << "No ribbon segments found for AcdId: " << id.id() << " with xOrientation "
                << xOrient << " Face: " << faceArr[iFace] << endreq;
	    return false;
        }
 
	// figure out segment is the "middle" one
	int toUse = (ribbonSegmentVolIds.size()-1) / 2;
	const idents::VolumeIdentifier& volId = ribbonSegmentVolIds[toUse];
	
	StatusCode sc;
	switch ( iFace ) {
	case 0:
	case 3:
	  // top
	  sc = m_glastDetSvc->getTransform3DByID(volId, &topTransform);
	  break;
	case 1:
	case 4:
	  // minus side
	  sc = m_glastDetSvc->getTransform3DByID(volId, &minusSideTransform);
	  break;
	case 2:
	case 5:
	  // plus side
	  sc = m_glastDetSvc->getTransform3DByID(volId, &plusSideTransform);
	  break;
	}
	if (sc.isFailure() ) {
	  MsgStream   log( msgSvc(), name() );
	  log << MSG::WARNING << "Failed to get ribbon trasnformation" << endreq;
	  return false;
	} 
    } // end face for


  return true;
}

double AcdGeometrySvc::ribbonHalfWidth() const 
{
  // FIXME (do this right)
  return 3.0;
}

/// Given an AcdId, provide the tile size, center and corners
bool AcdGeometrySvc::fillTileData(const idents::AcdId& id, int iVol,
				  int& face,
				  std::vector<double>& dim, 
				  HepPoint3D& center,
				  HepPoint3D* corner)
{
  idents::AcdId& ncid = const_cast<idents::AcdId&>(id);
  idents::VolumeIdentifier volId = ncid.volId(iVol==1);
  StatusCode sc = getDetectorDimensions(volId,dim,center);
  face = volId[1];
  if (sc.isFailure() ) {
    //log << MSG::WARNING << "Failed to get transformation" << endreq;
    return false;
  } 
  sc = getCorners(dim,center,corner);
  if (sc.isFailure() ) {
    //log << MSG::WARNING << "Failed to get transformation" << endreq;
    return false;
  } 
  return true;
}

/// Given an AcdId, provide transform to tile frame
bool AcdGeometrySvc::fillTileTransform(const idents::AcdId& id, int iVol,
				       HepTransform3D& transform)
{
  idents::AcdId& ncid = const_cast<idents::AcdId&>(id);
  idents::VolumeIdentifier volId = ncid.volId(iVol==1);
  StatusCode sc = m_glastDetSvc->getTransform3DByID(volId, &transform);
  if (sc.isFailure() ) {
    //log << MSG::WARNING << "Failed to get trasnformation" << endreq;
    return false;
  } 
  return true;
}

/// Given an AcdId, provide positions of screw holes in local frame
bool AcdGeometrySvc::fillScrewHoleData(const idents::AcdId& /* id */, 
				       std::vector< HepPoint3D >& /* screwHoles */) 
{
  // FIXME (do something)
  return true;
}

/// Given an AcdId, provide information about which volume edges are shared
bool AcdGeometrySvc::fillTileSharedEdgeData(const idents::AcdId& id, 
					    const std::vector<double>& dim1, const std::vector<double>& dim2,
					    int& sharedEdge1, int& sharedEdge2,
					    float& sharedWidth1, float& sharedWidth2)
{
  
  sharedEdge1 = sharedEdge2 = -1;
  sharedWidth1 = sharedWidth2 = 0.;
  switch (id.id()) {
  case 0: case 1: case 2: case 3: case 4:
    sharedEdge1 = 3;
    sharedEdge2 = 1;
    sharedWidth1 = dim2[2];
    sharedWidth2 = dim1[1];
    break;
  case 40: case 41: case 42: case 43: case 44:
    sharedEdge1 = 1;
    sharedEdge2 = 1;
    sharedWidth1 = dim2[2];
    sharedWidth2 = dim1[1];
    break;
  default:
    ;
  }
  return true;
}



#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "src/AcdGeometrySvc.h"
#include "AcdUtil/AcdFrameUtil.h"

#include "CLHEP/Geometry/Transform3D.h"

#include "idents/TowerId.h"

#include <iostream>
#include <algorithm>

static const SvcFactory<AcdGeometrySvc> s_factory;
const ISvcFactory& AcdGeometrySvcFactory = s_factory;

namespace {   // local utilities
  // Given a NamedId (std::vector<std::pair<std::string>, unsigned) look
  // for match with string field, store corresponding unsigned value.
  // in arg.  Return false if not found
  bool findFieldVal(const IGlastDetSvc::NamedId& nid, 
                    const std::string& fieldName, unsigned& val) {
    IGlastDetSvc::NamedId::const_iterator  it = nid.begin();
    while (it != nid.end()) {
      if (it->first == fieldName) {
        val = it->second;
        return true;
      }
      ++it;
    }
    return false;
  }
}

//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

AcdGeometrySvc::AcdGeometrySvc(const std::string& name, 
                               ISvcLocator* pSvcLocator) 
                              : Service(name, pSvcLocator)
{   

    AcdFrameUtil::buildRotations();
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
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to get AcdGeometrySvc constants" << endreq;
        return sc;
    }

    sc = getDetectorListFromGeometry();
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to get AcdGeometrySvc detector list" << endreq;
        return sc;
    }
 
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
    int eLATACD, eACDTile, eACDRibbon, eMeasureX, eMeasureY;
    int eACDTopFace, eACDXNegFace, eACDYNegFace, eACDXPosFace, eACDYPosFace;
    double ribbonWidth;
    if (m_glastDetSvc->getNumericConstByName("xNum", 
                                             &m_numXtowers).isSuccess() &&
        m_glastDetSvc->getNumericConstByName("yNum", 
                                             &m_numYtowers).isSuccess() &&
        m_glastDetSvc->getNumericConstByName("eLATACD", 
                                             &eLATACD).isSuccess()  &&
        m_glastDetSvc->getNumericConstByName("eACDTile", 
                                             &eACDTile).isSuccess() &&
        m_glastDetSvc->getNumericConstByName("eACDRibbon", 
                                             &eACDRibbon).isSuccess() &&
        m_glastDetSvc->getNumericConstByName("eMeasureX", 
                                             &eMeasureX).isSuccess() &&
        m_glastDetSvc->getNumericConstByName("eMeasureY", 
                                             &eMeasureY).isSuccess() &&
        m_glastDetSvc->getNumericConstByName("eACDTopFace", 
                                             &eACDTopFace).isSuccess() &&
        m_glastDetSvc->getNumericConstByName("eACDXNegFace", 
                                             &eACDXNegFace).isSuccess() &&
        m_glastDetSvc->getNumericConstByName("eACDYNegFace", 
                                             &eACDYNegFace).isSuccess() &&
        m_glastDetSvc->getNumericConstByName("eACDXPosFace", 
                                             &eACDXPosFace).isSuccess() &&
        m_glastDetSvc->getNumericConstByName("eACDYPosFace", 
                                             &eACDYPosFace).isSuccess() &&
        m_glastDetSvc->getNumericConstByName("ribbonWidth", 
                                             &m_ribbonHalfWidth).isSuccess() )
     {
       sc = StatusCode::SUCCESS;
       m_eLATACD = (unsigned) eLATACD;
       m_eACDTile = (unsigned) eACDTile;
       m_eACDRibbon = (unsigned) eACDRibbon;
       m_eMeasureX = (unsigned) eMeasureX;
       m_eMeasureY = (unsigned) eMeasureY;
       m_eACDTopFace = (unsigned) eACDTopFace;
       m_eACDXNegFace = (unsigned) eACDXNegFace;
       m_eACDYNegFace = (unsigned) eACDYNegFace;
       m_eACDXPosFace = (unsigned) eACDXPosFace;
       m_eACDYPosFace = (unsigned) eACDYPosFace;
       m_ribbonHalfWidth = ribbonWidth / 2;
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
        MsgStream log(msgSvc(), name());
        log << MSG::WARNING << "Failed to retrieve Shape by Id: " 
            << volId.name() << endreq;
        return sc;
    }
    HepGeom::Transform3D transform;
    sc = m_glastDetSvc->getTransform3DByID(volId, &transform);

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

    // Now get the TileDims
    const AcdTileDim* tile100 = m_geomMap.getTile(id100,*this);
    const AcdTileDim* tile104 = m_geomMap.getTile(id104,*this);
    const AcdTileDim* tile200 = m_geomMap.getTile(id200,*this);
    const AcdTileDim* tile204 = m_geomMap.getTile(id204,*this);
    const AcdTileDim* tile300 = m_geomMap.getTile(id300,*this);
    const AcdTileDim* tile304 = m_geomMap.getTile(id304,*this);
    const AcdTileDim* tile400 = m_geomMap.getTile(id400,*this);
    const AcdTileDim* tile404 = m_geomMap.getTile(id404,*this);
    const AcdTileDim* tile040 = m_geomMap.getTile(id040,*this);
    const AcdTileDim* tile430 = m_geomMap.getTile(id430,*this);

    // Check to see if these volIds exist in the geometry
    if ( tile100 == 0 )  sc = StatusCode::FAILURE;
    if ( tile104 == 0 ) sc = StatusCode::FAILURE;
    if ( tile200 == 0 ) sc = StatusCode::FAILURE;
    if ( tile204 == 0 ) sc = StatusCode::FAILURE;
    if ( tile300 == 0 ) sc = StatusCode::FAILURE;
    if ( tile304 == 0 ) sc = StatusCode::FAILURE;
    if ( tile400 == 0 ) sc = StatusCode::FAILURE;
    if ( tile404 == 0 ) sc = StatusCode::FAILURE;
    if ( tile040 == 0 ) sc = StatusCode::FAILURE;
    if ( tile430 == 0 ) sc = StatusCode::FAILURE;

    // If we fail to find any of the tiles, we skip calculation of the corner
    // gap variable
    if (sc.isFailure()) {
        log << MSG::WARNING << "Non-flight ACD Geometry will not calculate"
            << " corner gap rays" << endreq;
	// this is ok if we have the BeamTest release
	return StatusCode::SUCCESS;
    }

    // Determine the extent of all corner gap rays in Z, using the a top corner
    // tile and a bottom side row tile.
    double zMax = tile040->tileCenter().z();
    double zMin = tile430->corner()[0].z();    

    // Always use the top corner.  In local frame this is either -x,+y[1] or +x,+y[2]

    // Construct Ray for 1st Corner Gap, -x,-y corner
    // Use corners:  tile 100 +x,  tile 200 -x
    HepPoint3D corner0ref;

    AcdFrameUtil::getMidpoint(tile100->corner()[2],tile200->corner()[1],corner0ref);

    m_cornerGapStartPoint[0] = Point(corner0ref.x(),corner0ref.y(),zMin);
    m_cornerGapEndPoint[0] = Point(corner0ref.x(),corner0ref.y(),zMax);
    m_cornerGapVec[0] = Vector(m_cornerGapEndPoint[0] - m_cornerGapStartPoint[0]);

    log << MSG::DEBUG << "Corner1Ray:  Pos:(" << m_cornerGapStartPoint[0].x() 
        << ", " << m_cornerGapStartPoint[0].y() << ", " 
        << m_cornerGapStartPoint[0].z() << ") Vec:( "
        << m_cornerGapVec[0].x() << ", " << m_cornerGapVec[0].y() << ", " 
        << m_cornerGapVec[0].z() << ") mag: " << m_cornerGapVec[0].mag() << endreq;

    // Construct Ray for 2nd Corner Gap, -x, +y corner
    // Use corners:  tile 204 +x,  tile 300 -x
    HepPoint3D corner1ref;
    AcdFrameUtil::getMidpoint(tile204->corner()[2],tile300->corner()[1],corner1ref);

    m_cornerGapStartPoint[1] = Point(corner1ref.x(),corner1ref.y(),zMin);
    m_cornerGapEndPoint[1] = Point(corner1ref.x(),corner1ref.y(),zMax);
    m_cornerGapVec[1] = Vector(m_cornerGapEndPoint[1] - m_cornerGapStartPoint[1]);
    log << MSG::DEBUG << "Corner2Ray:  Pos:(" << m_cornerGapStartPoint[1].x() 
        << ", " << m_cornerGapStartPoint[1].y() << ", " 
        << m_cornerGapStartPoint[1].z() << ") Vec:( "
        << m_cornerGapVec[1].x() << ", " << m_cornerGapVec[1].y() << ", " 
        << m_cornerGapVec[1].z() << ") mag: " 
        << m_cornerGapVec[1].mag() << endreq;


    // Construct Ray for 3rd Corner Gap. +x, +y
    // Use corners:  tile 304 +x,  tile 404 -x
    HepPoint3D corner2ref;
    AcdFrameUtil::getMidpoint(tile304->corner()[2],tile404->corner()[1],corner2ref);

    m_cornerGapStartPoint[2] = Point(corner2ref.x(),corner2ref.y(),zMin);
    m_cornerGapEndPoint[2] = Point(corner2ref.x(),corner2ref.y(),zMax);
    m_cornerGapVec[2] = Vector(m_cornerGapEndPoint[2] - m_cornerGapStartPoint[2]);
    log << MSG::DEBUG << "Corner3Ray:  Pos:(" << m_cornerGapStartPoint[2].x() 
        << ", " << m_cornerGapStartPoint[2].y() << ", " 
        << m_cornerGapStartPoint[2].z() << ") Vec:( "
        << m_cornerGapVec[2].x() << ", " << m_cornerGapVec[2].y() << ", " 
        << m_cornerGapVec[2].z() << ") mag: " 
        << m_cornerGapVec[2].mag() << endreq;


    // Construct Ray for 4th Corner Gap, +x. -y
    // Use corners:  tile 400 -x,  tile 104 +x
    HepPoint3D corner3ref;
    AcdFrameUtil::getMidpoint(tile400->corner()[2],tile104->corner()[1],corner3ref);

    m_cornerGapStartPoint[3] = Point(corner3ref.x(),corner3ref.y(),zMin);
    m_cornerGapEndPoint[3] = Point(corner3ref.x(),corner3ref.y(),zMax);
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


  return retVal;
}

double AcdGeometrySvc::ribbonHalfWidth() const 
{
  // FIXME (do this right)
  return m_ribbonHalfWidth;
  //  return 3.0;
}

/// Given an AcdId, provide the tile size, center and corners
bool AcdGeometrySvc::fillTileData(const idents::AcdId& id, int iVol,
				  HepTransform3D& transformToLocal,
				  std::vector<double>& dim, 
				  HepPoint3D& center,
				  HepPoint3D* corner)
{

  MsgStream  log( msgSvc(), name() );

  // Build the VolumeIdentifier for this tile section
  idents::AcdId& ncid = const_cast<idents::AcdId&>(id);
  bool bent = iVol==1 ? true : false;
  const idents::VolumeIdentifier volId = ncid.volId(bent);

  // Get the reference frame enum
  AcdFrameUtil::AcdReferenceFrame frameId = getReferenceFrame(volId);
  if ( frameId == AcdFrameUtil::FRAME_NONE ) {
    log << MSG::WARNING << "Failed to retrieve Frame by Id: " 
	<< volId.name() << endreq;
    return false;
  }

  // dimensions in global frame
  std::vector<double> globalDim; 
  std::string str;

  // Now get the shape.
  // Note that this is expressed in the "GEANT" frame, 
  // which has only minimal rotations about X or Y axis for the side tiles
  StatusCode sc = m_glastDetSvc->getShapeByID(volId, &str, &globalDim);
  if ( sc.isFailure() ) {        
    log << MSG::WARNING << "Failed to retrieve Shape by Id: " 
	<< volId.name() << endreq;
    return false;
  } 
  
  // Make the dimension vector in the local frame;
  AcdFrameUtil::transformDimensionVector(frameId,globalDim,dim);

  // Get the transform from GEANT frame to the GLOBAL frame.  
  // Note that this includes the translation and is expressed in the global frame   
  HepGeom::Transform3D geantToGlobal;
  sc = m_glastDetSvc->getTransform3DByID(volId, &geantToGlobal);
  if (sc.isFailure() ) {
    log << MSG::WARNING << "Failed to get transformation: " 
	<< volId.name() << endreq;
    return false;
  }

  // Find the center of the tile, 
  HepPoint3D origin;
  center = geantToGlobal * origin;

  HepGeom::Transform3D globalToGeant = geantToGlobal.inverse();

  // Get the major rotations that flip axes around and all
  const HepTransform3D& rotationToLocal = AcdFrameUtil::getRotationToLocal(frameId);
  const HepTransform3D& rotationToGeant = AcdFrameUtil::getRotationToGeant(frameId);

  transformToLocal = rotationToLocal * globalToGeant;
  HepGeom::Transform3D transformToGlobal =  geantToGlobal * rotationToGeant;

  HepGeom::Transform3D check = transformToLocal * transformToGlobal;

  // Make the half-vectors (center to edge of tile)
  const HepVector3D xVectorLocal(dim[0]/2.,0.,0.);
  const HepVector3D yVectorLocal(0.,dim[1]/2.,0.);
  const HepVector3D xVectorGlobal = transformToGlobal* xVectorLocal;
  const HepVector3D yVectorGlobal = transformToGlobal* yVectorLocal;
  
  AcdFrameUtil::getCornersSquare(center,xVectorGlobal,yVectorGlobal,corner);

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

AcdFrameUtil::AcdReferenceFrame 
AcdGeometrySvc::getReferenceFrame(const idents::VolumeIdentifier &volId) {

    using idents::VolumeIdentifier;
    IGlastDetSvc::NamedId nid = m_glastDetSvc->getNamedId(volId);

    unsigned val;
    unsigned face;
    bool tile;

    if (!findFieldVal(nid, "fLATObjects", val)) return AcdFrameUtil::FRAME_NONE;
    if (val != m_eLATACD) return AcdFrameUtil::FRAME_NONE;

    if (!findFieldVal(nid, "fACDFace", face)) return AcdFrameUtil::FRAME_NONE;
    if (!findFieldVal(nid, "fACDCmp", val)) return AcdFrameUtil::FRAME_NONE;
    tile = (val == m_eACDTile);

    if (tile) {  // simple except for bent pieces
        if (face ==  m_eACDXNegFace) return AcdFrameUtil::FRAME_MINUSX;
        else if (face == m_eACDYNegFace) return AcdFrameUtil::FRAME_MINUSY;
        else if (face == m_eACDXPosFace) return AcdFrameUtil::FRAME_PLUSX;
        else if (face ==  m_eACDYPosFace) return AcdFrameUtil::FRAME_PLUSY;
        else if (face == m_eACDTopFace)            {
            if (!findFieldVal(nid, "fTileSeg", val)) return AcdFrameUtil::FRAME_NONE;
            if (val == 0) return AcdFrameUtil::FRAME_TOP;   // main tile piece
            if (!findFieldVal(nid, "fRow", val)) return AcdFrameUtil::FRAME_NONE;
            if (val == 0) return AcdFrameUtil::FRAME_MINUSY;
            else if (val == 4) return AcdFrameUtil::FRAME_PLUSY_YDWN;
            else return AcdFrameUtil::FRAME_NONE;
        }
        else return AcdFrameUtil::FRAME_NONE;
    }
    else if (val != m_eACDRibbon) return AcdFrameUtil::FRAME_NONE;
  
    // Ribbons...eek!
    if (!findFieldVal(nid, "fMeasure", val)) return AcdFrameUtil::FRAME_NONE;
    unsigned ribbon, segNum;
    if (!findFieldVal(nid, "fRibbon", ribbon)) return AcdFrameUtil::FRAME_NONE;
    if (!findFieldVal(nid, "fRibbonSegment", segNum)) return AcdFrameUtil::FRAME_NONE;

    bool increasing;
    std::vector<double> dims;
    HepPoint3D cm;
    if (getDetectorDimensions(volId, dims, cm).isFailure()) return AcdFrameUtil::FRAME_NONE;

    //    HepGeom::Transform3D transform;
    //    m_glastDetSvc->getTransform3DByID(volId, &transform);
    //    HepVector3D dimsTransformed(dims[0], dims[1], dims[2]);
    //    dimsTransformed = transform * dimsTransformed;
    // More generally might need to use transformed dimensions in 
    //   following, but for actual geometry in use it isn't necessary
    bool measY = true;

    if (val == m_eMeasureX) {  // width is X dimension.  
        measY = false;
        if (face == m_eACDTopFace) {
          //            if (dims[2] <= dims[1]) return FRAME_XMEAS;
            if (dims[2] <= dims[1]) return AcdFrameUtil::FRAME_XMEAS;
            increasing = true;
        }
        else if (face == m_eACDYNegFace) {
          //  if (dims[2] >= dims[1]) return FRAME_MINUSY;  // vert
          if (dims[2] >= dims[1]) return AcdFrameUtil::FRAME_MINUSY;  // it's vertical 
            increasing = true;
        }
        else if (face == m_eACDYPosFace) {
            if (dims[2] >= dims[1]) return AcdFrameUtil::FRAME_PLUSY_YDWN; // it's vertical
            increasing = false;
        }
        else return AcdFrameUtil::FRAME_NONE;
    }
    else if (val = m_eMeasureY) {  // width is Y dimension

        // Y-ribbons go straight across top
        if (face == m_eACDTopFace) return AcdFrameUtil::FRAME_YMEAS; 
        else if (face == m_eACDXNegFace) {
            if (dims[2] >= dims[0]) return AcdFrameUtil::FRAME_MINUSX;  // it's vertical 
            increasing = true;
        }
        //            break;
        else if (face == m_eACDXPosFace) {
            if (dims[2] >= dims[0]) return AcdFrameUtil::FRAME_PLUSX_YDWN; // it's vertical
            increasing = false;
        }
        //            break;
        else return AcdFrameUtil::FRAME_NONE;
    }
    // deal with short guys.
    std::vector<VolumeIdentifier> segs;

    m_glastDetSvc->orderRibbonSegments(segs, face, ribbon, measY, 
                                       increasing);
    VolIdIter ourIt = std::find(segs.begin(), segs.end(), volId);
    // If our seg wasn't found, give up
    if (ourIt == segs.end())  return AcdFrameUtil::FRAME_NONE;
    // Our seg really shouldn't be first or last
    if ((ourIt == segs.begin()) || ((ourIt + 1) == segs.end()) )
        return AcdFrameUtil::FRAME_NONE; 

    VolIdIter before = ourIt - 1;
    VolIdIter after = ourIt + 1;

    HepPoint3D xTbefore, xTafter;
    if (getDimensions(*before, dims, xTbefore).isFailure()) return AcdFrameUtil::FRAME_NONE;
    if (getDimensions(*after, dims, xTafter).isFailure()) return AcdFrameUtil::FRAME_NONE;

    // Depending on face, compare x, y or z dimensions of cm
    // to see if we're headed forwards or backward
    // Y-measuring, side
    if ((face ==  m_eACDXNegFace) || (face ==  m_eACDXPosFace)) {
        if (xTbefore.x() <= xTafter.x()) return AcdFrameUtil::FRAME_YMEAS;
        else return AcdFrameUtil::FRAME_YMEAS_ZROT180;
    }

    // X-measuring, side
    else if ((face ==  m_eACDYNegFace) || (face == m_eACDYPosFace)) {
        if (xTbefore.y() <= xTafter.y()) return AcdFrameUtil::FRAME_XMEAS;
        else return AcdFrameUtil::FRAME_XMEAS_ZROT180;
    }
    else {   // must be top
        if (xTbefore.z() <= xTafter.z()) return AcdFrameUtil::FRAME_MINUSY;
        else return AcdFrameUtil::FRAME_PLUSY_YDWN;
    }
}


#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "src/AcdGeometrySvc.h"
#include "AcdUtil/AcdFrameUtil.h"

#include "CLHEP/Geometry/Transform3D.h"

#include "xmlBase/Dom.h"
#include "xmlBase/XmlParser.h"
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMElement.hpp>

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
    declareProperty("AlignmentFile",    m_alignFileName = "");
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

    sc = readAlignmentFile();
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to read alignment" << endreq;
        return sc;
    }
     
    log << MSG::INFO << "AcdGeometrySvc successfully initialized" << endreq;
    return StatusCode::SUCCESS;

}


StatusCode AcdGeometrySvc::readAlignmentFile() {

  XERCES_CPP_NAMESPACE_USE;

  std::string alignFile = m_alignFileName;
  if ( alignFile.size() == 0 ) return StatusCode::SUCCESS;
  
  MsgStream log(msgSvc(), name());
  log << MSG::WARNING << "Reading ACD aligment from " << alignFile << endreq;

  
  xmlBase::XmlParser parser;
  DOMDocument* doc(0);
  std::vector<DOMElement*> tileElems;
  try {
    doc = parser.parse(alignFile.c_str());
    if ( doc == 0 ) return StatusCode::FAILURE;
    const DOMElement* root = doc->getDocumentElement();
    xmlBase::Dom::getDescendantsByTagName(root,"tile",tileElems);
    for ( std::vector<DOMElement*>::const_iterator itr = tileElems.begin();
	  itr != tileElems.end(); itr++ ) {
      int tileId = xmlBase::Dom::getIntAttribute(*itr,"tileId");
      std::vector<DOMElement*> volElems;
      xmlBase::Dom::getDescendantsByTagName(*itr,"volume",volElems);
      idents::AcdId acdId(tileId);
      for (  std::vector<DOMElement*>::const_iterator itr2 = volElems.begin();
	     itr2 != volElems.end(); itr2++ ) {
	int volId = xmlBase::Dom::getIntAttribute(*itr2,"iVol");
	const idents::VolumeIdentifier& vId = acdId.volId( volId > 0 );
	DOMElement* alignElem =  xmlBase::Dom::findFirstChildByName(*itr2,"acdAlign");
	double centerX = xmlBase::Dom::getDoubleAttribute(alignElem,"centerX");
	double centerY = xmlBase::Dom::getDoubleAttribute(alignElem,"centerY");
	double sizeX = xmlBase::Dom::getDoubleAttribute(alignElem,"sizeX");
	double sizeY = xmlBase::Dom::getDoubleAttribute(alignElem,"sizeY");
	AcdUtil::AcdVolumeAlignment align(centerX,centerY,sizeX,sizeY);
	m_geomMap.putAlignment(vId,align);
      }
    } 
  } catch ( xmlBase::DomException& e ) {
    log << MSG::ERROR << "Xml parsing of " << alignFile << " failed with exception: " << std::endl 
	<< "\t" << e.getMsg() << endreq;
    return StatusCode::FAILURE;
  }
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
                                             &ribbonWidth).isSuccess() )
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


bool AcdGeometrySvc::fillRibbonData(const idents::AcdId& id,
				    std::vector<Ray>& minusSideRays,
				    std::vector<Ray>& topRays,
				    std::vector<Ray>& plusSideRays, 
				    HepTransform3D& minusSideTransform,
				    HepTransform3D& topTransform,
				    HepTransform3D& plusTransform) {
  // Purpose and Method:  Fill the three supplied vector of Rays.  The Rays are constructed from the ribbon segments
  //    associated with AcdId id.  
  MsgStream   log( msgSvc(), name() );
  
  log << MSG::DEBUG << "Filling data from ribbon " << id.id() << endreq;

  minusSideRays.clear();
  topRays.clear();
  plusSideRays.clear();

  typedef enum {
    ribbonX = 5,
    ribbonY = 6
  } ribbonOrient;

  bool retVal = true;
  if (!id.ribbon()) return false;
  
  unsigned int ribbonNum = id.ribbonNum();
  unsigned int ribbonOrientation = id.ribbonOrientation();
  
  static const unsigned int nRibbonXSideSegsUsed(5);
  static const unsigned int nRibbonXTopSegsUsed(1);
  
  static const unsigned int nRibbonYSideSegsUsed(7);
  static const unsigned int nRibbonYTopSegsUsed(5);
  
  // These numberings are in terms of the geometrical order of the segments, 
  // not their segement field in the volume identifier  
  static const unsigned int ribbonXSideSegs[nRibbonXSideSegsUsed] = { 0,1,2,3,4 };     // use all the segments
  static const unsigned int ribbonXTopSegs[nRibbonXTopSegsUsed] = { 0 };               // only one segment

  static const unsigned int ribbonYSideSegs[nRibbonYSideSegsUsed] = { 0,1,2,3,4,5,6 }; // use all the segments
  static const unsigned int ribbonYTopSegs[nRibbonYTopSegsUsed] = { 0,2,4,6,8 };       // use only the long segments
  
  static const unsigned int xFaces[3] = { 1,0,3 }; // which faces of the detector the segment lie along
  static const unsigned int yFaces[3] = { 2,0,4 };

  static const unsigned int xRefSeg[3] = { 4,0,0 }; // reference segements for "ribbon frame"
  static const unsigned int yRefSeg[3] = { 4,2,2 }; 

  // orientation of ribbon 
  bool xOrient = ribbonOrientation == ribbonX ? true : false;

  static const Point nullPoint;
  static const Vector nullVector;
  static const Ray nullRay(nullPoint,nullVector);
 
  if ( xOrient ) {
    minusSideRays.resize(nRibbonXSideSegsUsed,nullRay);
    topRays.resize(nRibbonXTopSegsUsed,nullRay);
    plusSideRays.resize(nRibbonXSideSegsUsed,nullRay);
  } else {
    minusSideRays.resize(nRibbonYSideSegsUsed,nullRay);
    topRays.resize(nRibbonYTopSegsUsed,nullRay);
    plusSideRays.resize(nRibbonYSideSegsUsed,nullRay);
  }

  // Loops over set of faces
  for (unsigned iFace = 0; iFace < 3; iFace++) {
      
    // Get the ribbon segments from the glastDetSvc
    std::vector<idents::VolumeIdentifier> ribbonSegmentVolIds;
    m_glastDetSvc->orderRibbonSegments(ribbonSegmentVolIds,
				       xOrient ? xFaces[iFace] : yFaces[iFace], ribbonNum, xOrient, true);

    
    // OK, now loop over the relevent segments, first we need to figure out which they are
    const unsigned int* segmentIndex = xOrient ? 
      ( iFace == 1 ? ribbonXTopSegs : ribbonXSideSegs ) :
      ( iFace == 1 ? ribbonYTopSegs : ribbonYSideSegs );
    unsigned int nSegment = xOrient ? 
      ( iFace == 1 ? nRibbonXTopSegsUsed : nRibbonXSideSegsUsed ) :
      ( iFace == 1 ? nRibbonYTopSegsUsed : nRibbonYSideSegsUsed );

    unsigned int checkRefSegment = xOrient ? xRefSeg[iFace] : yRefSeg[iFace];

    for ( unsigned int iSeg(0); iSeg < nSegment; iSeg++ ) {

      const idents::VolumeIdentifier& volId = ribbonSegmentVolIds[  segmentIndex[iSeg]  ];
            
      // Make the dimension vector in the local frame;
      std::vector<double> dim;
      HepGeom::Transform3D transformToLocal;
      HepPoint3D center;
      HepVector3D xVectorGlobal;
      HepVector3D xVectorDummy;
      HepVector3D yVectorGlobal;
      
      StatusCode sc = getTransformAndLocalVectors(volId,dim,transformToLocal,center,xVectorGlobal,xVectorDummy, yVectorGlobal);

      if ( sc.isFailure() ) {        
	log << MSG::ERROR << "Failed to handle transformations for ribbon volume: " 
	    << volId.name() << endreq;
	return sc;
      } 

      HepPoint3D start = center - yVectorGlobal;
      HepPoint3D end = center + yVectorGlobal;
      HepVector3D vect = 2 * yVectorGlobal;

      Point startPoint(start.x(),start.y(),start.z());
      Vector rayVector(vect.x(),vect.y(),vect.z());

      unsigned rayIndex(0);
      switch (iFace) {
      case 0: 
	rayIndex = iSeg;
	minusSideRays[rayIndex] = Ray(startPoint,rayVector);
	minusSideRays[rayIndex].setArcLength(rayVector.mag());
	break;
      case 1: 
	rayIndex = iSeg;
	topRays[rayIndex] = Ray(startPoint,rayVector);
	topRays[rayIndex].setArcLength(rayVector.mag());
	break;
      case 2: 
	rayIndex = nSegment - (iSeg+1);
	plusSideRays[rayIndex] = Ray(startPoint,rayVector);
	plusSideRays[rayIndex].setArcLength(rayVector.mag());
	break;
      default:
	return false;
      }
     
      if ( iSeg == checkRefSegment ) {
	switch (iFace) {
	case 0:
	  minusSideTransform = transformToLocal;
	  break;
	case 1:
	  topTransform = transformToLocal;
	  break;
	case 2:
	  plusTransform = transformToLocal;
	  break;
	default:
	  return false;	  
	}
      }
      
      log << MSG::DEBUG << volId.name() << ' ' << (xOrient ? xFaces[iFace] : yFaces[iFace]) << ' ' << iSeg << ' ' << (xOrient ? 'X' : 'Y') 
	  << ' ' << ribbonNum << ' ' << (iSeg == checkRefSegment ? "REF" : "") << std::endl
	  << "cen: " << center << ", " << std::endl
	  << "dim: " << dim[0] << ", " << dim[1] << ", " << dim[2] << std::endl
	  << "startPos: (" << startPoint.x() << ", " << startPoint.y() << ", " << startPoint.z() << ")" 
	  << " endPos: ( " << end.x() << ", " << end.y() << ", " << end.z() << ")" << endreq;
      

    }

  }

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
  
  log << MSG::DEBUG << "Filling data from tile " << id.id() << ':' << iVol << endreq;

  // Build the VolumeIdentifier for this tile section
  idents::AcdId& ncid = const_cast<idents::AcdId&>(id);
  bool bent = iVol==1 ? true : false;
  const idents::VolumeIdentifier volId = ncid.volId(bent);

  HepVector3D x1VectorGlobal;
  HepVector3D x2VectorGlobal;
  HepVector3D yVectorGlobal;
  
  StatusCode sc = getTransformAndLocalVectors(volId,dim,transformToLocal,center,
					      x1VectorGlobal,x2VectorGlobal,yVectorGlobal);  

  if ( sc.isFailure() ) {        
    log << MSG::ERROR << "Failed to handle transformations for tile volume: " 
      	<< volId.name() << endreq;
    return sc;
  } 

  double x2(0.), xdelta(0.);
  if ( dim.size() == 3 ) {   
    AcdFrameUtil::getCornersSquare(center,x1VectorGlobal,yVectorGlobal,corner);
  } else {
    AcdFrameUtil::getCornersTrap(center,x1VectorGlobal,x2VectorGlobal,yVectorGlobal,corner);
    // Also grab a couple of things for the printout
    x2 = dim[3];
    xdelta = dim[4];
  }
  
  
  log << MSG::DEBUG << volId.name() << ' ' 
      << "cen: " << center << ", " << std::endl
      << "dim: " << dim[0] << ", " << dim[1] << ", " << dim[2] << ", " << x2 << ", " << xdelta << std::endl
      << "x1v: " << x1VectorGlobal << std::endl
      << "x2v: " << x2VectorGlobal << std::endl
      << "yv: "  << yVectorGlobal << std::endl
      << "corner: (" << corner[0] << ',' << corner[1] << ',' << corner[2] << ',' <<  corner[3] << ',' << ")" << endreq;      

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
  
  // edges are in order (-x, +y, +x, -y)
  sharedEdge1 = sharedEdge2 = -1;
  sharedWidth1 = sharedWidth2 = 0.;
  switch (id.id()) {
  case 0: case 1: case 2: case 3: case 4:
    // top -Y side tiles
    // bent piece orientation is (+x,+z,-y)
    sharedEdge1 = 3;   // top piece shares local -y edge (-y)
    sharedEdge2 = 1;   // bent piece shares local +y edge (+z)
    sharedWidth1 = dim2[1];  // extra width is local y  
    sharedWidth2 = dim1[1];  // extra width is local y
    break;
  case 40: case 41: case 42: case 43: case 44:
    // top +Y side tiles
    // bent piece orientation is (+x,-z,+y)
    sharedEdge1 = 1;   // top piece shares local +y edge (+y)
    sharedEdge2 = 3;   // bent pieces shares local -y edge(+z)
    sharedWidth1 = dim2[1];  // extra width is local y
    sharedWidth2 = dim1[1];  // extra width is local y
    break; 
 default:
    ;
  }
  return true;
}

 AcdFrameUtil::AcdReferenceFrame
AcdGeometrySvc::getReferenceFrame(const idents::VolumeIdentifier &volId) const {

    using idents::VolumeIdentifier;
    IGlastDetSvc::NamedId nid = m_glastDetSvc->getNamedId(volId);

    unsigned val;
    unsigned face;
    bool tile;

    // Make sure it is part of the ACD
    if (!findFieldVal(nid, "fLATObjects", val)) return AcdFrameUtil::FRAME_NONE;
    if (val != m_eLATACD) return AcdFrameUtil::FRAME_NONE;

    // Check the face and make if it is a tile or ribbon
    if (!findFieldVal(nid, "fACDFace", face)) return AcdFrameUtil::FRAME_NONE;
    if (!findFieldVal(nid, "fACDCmp", val)) return AcdFrameUtil::FRAME_NONE;
    tile = (val == m_eACDTile);

    if (tile) {  // simple except for bent pieces
      // Everything but the top just depends on the face
      // if (face ==  m_eACDXNegFace) return AcdFrameUtil::FRAME_MINUSX;
      // else if (face == m_eACDYNegFace) return AcdFrameUtil::FRAME_MINUSY;
      // else if (face == m_eACDXPosFace) return AcdFrameUtil::FRAME_PLUSX;
      // else if (face ==  m_eACDYPosFace) return AcdFrameUtil::FRAME_PLUSY;
      if ((face == m_eACDXNegFace) || (face == m_eACDYNegFace) || 
          (face == m_eACDXPosFace) || (face == m_eACDYPosFace)) 
        return AcdFrameUtil::FRAME_TOP;
      // For the to check to see if it is a bent piece
      else if (face == m_eACDTopFace) {
	if (!findFieldVal(nid, "fTileSeg", val)) return AcdFrameUtil::FRAME_NONE;
	// Is main tile piece, return FRAME_TOP
	if (val == 0) return AcdFrameUtil::FRAME_TOP;   
	// It is the bent piece, return the frame that is an extension of the rest of the tile
	if (!findFieldVal(nid, "fRow", val)) return AcdFrameUtil::FRAME_NONE;
	// Row 0 has Y going up the side
	if (val == 0) return AcdFrameUtil::FRAME_MINUSY;
	// Row 4 has Y going down the side
	else if (val == 4) return AcdFrameUtil::FRAME_PLUSY_YDWN;
	else return AcdFrameUtil::FRAME_NONE;
      }
      else return AcdFrameUtil::FRAME_NONE;
    }
    else if (val != m_eACDRibbon) return AcdFrameUtil::FRAME_NONE;
  
    // Ribbons...eek!
    if (!findFieldVal(nid, "fMeasure", val)) return AcdFrameUtil::FRAME_NONE;
    unsigned segNum;
    if (!findFieldVal(nid, "fRibbonSegment", segNum)) return AcdFrameUtil::FRAME_NONE;

    static const AcdFrameUtil::AcdReferenceFrame 
      xTopFrames[12] = { AcdFrameUtil::FRAME_NONE,
			 AcdFrameUtil::FRAME_XMEAS,AcdFrameUtil::FRAME_XMEAS,AcdFrameUtil::FRAME_XMEAS,
			 AcdFrameUtil::FRAME_XMEAS,AcdFrameUtil::FRAME_XMEAS,
			 AcdFrameUtil::FRAME_NONE,AcdFrameUtil::FRAME_NONE,
			 AcdFrameUtil::FRAME_PLUSY_YDWN,AcdFrameUtil::FRAME_PLUSY_YDWN, 
			 AcdFrameUtil::FRAME_MINUSY,AcdFrameUtil::FRAME_MINUSY };

    if ( val == m_eMeasureX ) {
      if ( face == m_eACDTopFace ) return xTopFrames[segNum];
      else if ( face == m_eACDYNegFace ) return AcdFrameUtil::FRAME_MINUSY;
      else if ( face == m_eACDYPosFace ) return AcdFrameUtil::FRAME_PLUSY_YDWN;
      else return AcdFrameUtil::FRAME_NONE;
    } else if ( val == m_eMeasureY ) {
      if ( face == m_eACDTopFace) return AcdFrameUtil::FRAME_YMEAS; 
      else if (face == m_eACDXNegFace) return AcdFrameUtil::FRAME_MINUSX;  
      else if (face == m_eACDXPosFace) return AcdFrameUtil::FRAME_PLUSX_YDWN; 
      else return AcdFrameUtil::FRAME_NONE;
    }  
    return AcdFrameUtil::FRAME_NONE;

}


StatusCode AcdGeometrySvc::getTransformAndLocalVectors(const idents::VolumeIdentifier &volId,
						       std::vector<double>& dim,
						       HepGeom::Transform3D& transformToLocal,
						       HepPoint3D& center,
						       HepVector3D& x1VectorGlobal,
						       HepVector3D& x2VectorGlobal,
						       HepVector3D& yVectorGlobal) const {

  MsgStream  log( msgSvc(), name() );

  // dimensions in global frame
  std::vector<double> globalDim; 
  std::string str;

  // Now get the shape.
  // Note that this is expressed in the "GEANT" frame, 
  // which has only minimal rotations about X or Y axis for the side tiles
  StatusCode sc = m_glastDetSvc->getShapeByID(volId, &str, &globalDim);

  if ( sc.isFailure() ) {        
    log << MSG::ERROR << "Failed to retrieve shape by Id: " 
	<< volId.name() << endreq;
    return sc;
  } 

  // Get the reference frame enum
  AcdFrameUtil::AcdReferenceFrame frameId = getReferenceFrame(volId);
  if ( frameId == AcdFrameUtil::FRAME_NONE ) {
    log << MSG::ERROR << "Failed to retrieve Frame by Id: " 
	<< volId.name() << endreq;
    return StatusCode::FAILURE;
  }

  // Get the alignment for this volume
  const AcdUtil::AcdVolumeAlignment& align = m_geomMap.getAlignment(volId);

  // Make the dimension vector in the local frame;
  if ( globalDim.size() == 3 ) {
    AcdFrameUtil::transformDimensionVector(frameId,globalDim,dim);
    if ( align.sizeX() > 0 ) {
      dim[0] = align.sizeX();
    }
    if ( align.sizeY() > 0 ) {
      dim[1] = align.sizeY();
    }    
  } else {    
    AcdFrameUtil::transformDimensionVectorTrap(frameId,globalDim,dim);
    if ( align.sizeX() > 0 ) {
      // FIXME: too complicated, drop it for now.
    }    
    if ( align.sizeY() > 0 ) {
      dim[1] = align.sizeY();
    }    
  }

  // Get the transform from GEANT frame to the GLOBAL frame.  
  // Note that this includes the translation and is expressed in the global frame   
  HepGeom::Transform3D geantToGlobal;
  sc = m_glastDetSvc->getTransform3DByID(volId, &geantToGlobal);
  if (sc.isFailure() ) {
    log << MSG::ERROR << "Failed to get transformation: " 
	<< volId.name() << endreq;
    return sc;
  }


  // Find the center of the volume, 
  HepPoint3D origin(align.centerX(),align.centerY(),0.);
  center = geantToGlobal * origin;


  // Build up the transformation to the local frame
  // Get the minor rotations and offsets that take the global frame to the geant reference
  HepGeom::Transform3D globalToGeant = geantToGlobal.inverse();
  // Get the major rotations that flip axes around and all
  const HepTransform3D& rotationToLocal = AcdFrameUtil::getRotationToLocal(frameId);  
  transformToLocal = rotationToLocal * globalToGeant;

  // Build up the transformation to the global frame
  // Get the major rotations that flip axes around and all
  const HepTransform3D& rotationToGeant = AcdFrameUtil::getRotationToGeant(frameId);
  HepGeom::Transform3D transformToGlobal =  geantToGlobal * rotationToGeant;

  // This is just here as a sanity check should be identity
  //HepGeom::Transform3D check = transformToLocal * transformToGlobal;

  // Make the half-vectors (center to edge of volume)  
  double xHalfLength = dim[0]/2;
  double yHalfLength = dim[1]/2;
  double xOffsetY = dim.size() == 3 ? 0. : dim[4]/2.;
  
  const HepVector3D x1VectorLocal(xHalfLength,0.,0.);
  x1VectorGlobal = transformToGlobal* x1VectorLocal;

  const HepVector3D yVectorLocal(xOffsetY,yHalfLength,0.);
  yVectorGlobal = transformToGlobal* yVectorLocal;

  if ( dim.size() > 3 ) {
    const HepVector3D x2VectorLocal(dim[3]/2,0.,0.);
    x2VectorGlobal = transformToGlobal* x2VectorLocal;
  }

  // done, return
  return sc;
  
}
			    

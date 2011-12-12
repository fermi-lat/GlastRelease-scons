#define TestAcdUtil_CXX

// File and Version Information
// $Header$
// Description:
// Test for AcdUtil class. 

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "GaudiKernel/ObjectVector.h"

#include "GaudiKernel/Algorithm.h"

#include "AcdUtil/IAcdGeometrySvc.h"
#include "AcdUtil/AcdDetectorList.h"
#include "AcdUtil/AcdRibbonDim.h"
#include "AcdUtil/AcdTileDim.h"
#include "AcdUtil/AcdFrameUtil.h"
#include "AcdUtil/IAcdFailureModeSvc.h"

#include "idents/AcdId.h"


/** @class TestAcdUtil
 * @brief AcdUtil test algorithm
 *
 * Exercise all of AcdUtil to be sure that the methods function properly.
 *
 * $Header$
 */

class TestAcdUtil : public Algorithm {

public:
  TestAcdUtil(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  MsgStream *m_log;

  IGlastDetSvc    *m_glastDetSvc;
  IAcdGeometrySvc *m_acdGeoSvc;
  IAcdFailureModeSvc *m_acdFailureSvc;

  void writeTileFrame(unsigned face, unsigned row, unsigned col, 
                      bool bent = false);
  void writeRibbonFrame(unsigned face, unsigned orient, unsigned rib, 
                        unsigned seg);
};

//static const AlgFactory<TestAcdUtil>  Factory;
//const IAlgFactory& TestAcdUtilFactory = Factory;
DECLARE_ALGORITHM_FACTORY(TestAcdUtil);

TestAcdUtil::TestAcdUtil(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator), m_log(0) {
    
}

StatusCode TestAcdUtil::initialize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    StatusCode sc = StatusCode::SUCCESS;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    m_glastDetSvc = 0;
    sc = service("GlastDetSvc", m_glastDetSvc, true);
    if (sc.isSuccess() ) {
        sc = m_glastDetSvc->queryInterface(IID_IGlastDetSvc, (void**)&m_glastDetSvc);
    }

    sc = service("AcdGeometrySvc", m_acdGeoSvc, true);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to find ACD Geometry service" << endreq;
        return sc;
    }

    sc = service("AcdFailureModeSvc", m_acdFailureSvc, true);
    if (sc.isFailure() ){
        log << MSG::ERROR << " Unable to find ACD Failure Mode service" << endreq;
        return sc;
    }

    // Testing AcdFailureModeSvc
    idents::AcdId id500(0,5,0,0);
    idents::AcdId id600(0,6,0,0);
    idents::AcdId id310(0,3,1,0);
    if (m_acdFailureSvc->matchAcdId(id500))
        log << MSG::INFO << "AcdId 500 found in AcdFailureModeSvc" << endreq;
    if (m_acdFailureSvc->matchAcdId(id600))
        log << MSG::INFO << "AcdId 600 found in AcdFailureModeSvc" << endreq;
    if (!m_acdFailureSvc->matchAcdId(id310)) 
        log << MSG::INFO << "AcdId 310 not found in AcdFailureModeSvc" << endreq;


    //m_acdGeoSvc->findCornerGaps();

    return sc;
}


StatusCode TestAcdUtil::execute() {
    // Purpose and Method:  Exercise the methods contained in AcdDigiUtil
    //  When sampling the poisson and gaussian distributions, store the results to
    //  the appropriate vector.

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    if (m_log == 0) m_log = new MsgStream(msgSvc(), name());

    log << MSG::INFO <<  "Number of Tiles: " << m_acdGeoSvc->numTiles()
        << " Number of Ribbons: " << m_acdGeoSvc->numRibbons() << endreq;

    const AcdUtil::AcdDetectorList& acdList =  m_acdGeoSvc->getDetectorList(); 

    AcdUtil::AcdDetectorList::const_iterator acdIt;

    unsigned int acdDetectorCounter = 0;
    for (acdIt = acdList.begin(); acdIt != acdList.end(); acdIt++) {
        ++acdDetectorCounter;
        idents::VolumeIdentifier volid = *acdIt;
        bool found = m_acdGeoSvc->findDetector(volid);
        if (!found) {
            log << MSG::WARNING << "Failed to find " << volid.name() << endreq;
            return StatusCode::FAILURE;
        }
	
	idents::AcdId acdId(volid);

	if ( acdId.tile() ) { 
	  // Test creation of AcdTileDim
	  AcdTileDim aTile(acdId,*m_acdGeoSvc);
	  int section = volid[5];

	  log << MSG::INFO << "AcdId: " << acdId.id() << " Dimensions: " 
	      << aTile.getSection(0)->m_dim[0] << ", " << aTile.getSection(0)->m_dim[1] << ", " << aTile.getSection(0)->m_dim[2] << endreq;

	  unsigned int iCorner;
	  for (iCorner = 0; iCorner<4; iCorner++) {
            log << "Corner " << iCorner << " " << aTile.getSection(0)->m_corners[iCorner].x() 
                << ", " << aTile.getSection(0)->m_corners[iCorner].y() << ", " << aTile.getSection(0)->m_corners[iCorner].z()
                << endreq;
	  }	  
	} else if (acdId.ribbon()) {
            // Test Creation of RibbonDim
	    AcdRibbonDim ribbonDim(acdId,*m_acdGeoSvc);
            //const HepPoint3D* start = ribbonDim.ribbonStart();
            //const HepPoint3D* end = ribbonDim.ribbonEnd();
            //const double* halfWidth = ribbonDim.halfWidth();
     
            log << MSG::DEBUG << "ID: " << acdId.id() << std::endl;
	      //<< " start:(" 
              //  << start[1].x() << "," << start[1].y() << "," << start[1].z()
              //  << ") end:(" << end[1].x() << "," << end[1].y() << ","
              //  << end[1].z() << ") halfWid: " << halfWidth[1] << std::endl;

	    std::vector<AcdRibbonSegment*> segs;
	    int topIdx(0), plusIdx(0);
            m_acdGeoSvc->fillRibbonData(acdId, segs, topIdx, plusIdx);
            log << MSG::DEBUG << "Trays: " << endreq;
            for ( std::vector<AcdRibbonSegment*>::const_iterator segIt = segs.begin(); segIt != segs.end(); segIt++) {
  	        const HepPoint3D& pos = (*segIt)->m_start;
	        const HepVector3D& dir = (*segIt)->m_vect;
                log << MSG::DEBUG << "Pos: ( " << pos.x() << ", " << pos.y() << ", " << pos.z() << ")" << endreq;
                log << MSG::DEBUG << "Dir: ( " << dir.x() << ", " << dir.y() << ", " << dir.z() << ")" << endreq;
            }
        }

    }


    log << MSG::INFO << "Number of Detectors:  " << acdDetectorCounter << endreq;
    // Test getReferenceFrame for a few different cases
    idents::VolumeIdentifier vid;
    idents::AcdId aid(0,0,0,2);
    vid = aid.volId();

    writeTileFrame(0, 0, 0);
    writeTileFrame(0, 0, 0, true);

    writeTileFrame(2, 1, 2);
    writeTileFrame(0, 2, 2);
    writeTileFrame(0, 4, 2);
    writeTileFrame(0, 4, 2, true);

    writeTileFrame(1, 4, 3);
    writeTileFrame(4, 0, 3);

    // Now try some ribbons
    log << "Y-measuring ribbons: " << endreq;
    writeRibbonFrame(0, 1, 2, 1);

    writeRibbonFrame(1, 1, 2, 1);
    writeRibbonFrame(1, 1, 3, 2);
    writeRibbonFrame(1, 1, 2, 3);
    writeRibbonFrame(1, 1, 3, 4);
    writeRibbonFrame(1, 1, 2, 6);
    writeRibbonFrame(1, 1, 3, 7);
    writeRibbonFrame(1, 1, 0, 8);

    writeRibbonFrame(3, 1, 3, 4);
    writeRibbonFrame(3, 1, 2, 6);
    writeRibbonFrame(3, 1, 3, 7);
    writeRibbonFrame(3, 1, 0, 8);

    log << "X-measuring top ribbons: " << endreq;
    writeRibbonFrame(0, 0, 2, 1);
    writeRibbonFrame(0, 0, 2, 2);
    writeRibbonFrame(0, 0, 2, 3);
    writeRibbonFrame(0, 0, 3, 4);
    writeRibbonFrame(0, 0, 1, 8);
    writeRibbonFrame(0, 0, 2, 9);
    writeRibbonFrame(0, 0, 3, 10);

    log << "X-measuring side ribbons: " << endreq;
    writeRibbonFrame(2, 0, 1, 1);
    writeRibbonFrame(2, 0, 1, 2);
    writeRibbonFrame(2, 0, 1, 3);
    writeRibbonFrame(2, 0, 1, 6);
    writeRibbonFrame(2, 0, 1, 7);
    writeRibbonFrame(2, 0, 1, 8);

    writeRibbonFrame(4, 0, 1, 1);
    writeRibbonFrame(4, 0, 1, 2);
    writeRibbonFrame(4, 0, 1, 3);
    writeRibbonFrame(4, 0, 1, 6);
    writeRibbonFrame(4, 0, 1, 7);
    writeRibbonFrame(4, 0, 1, 8);


    return sc;
}


StatusCode TestAcdUtil::finalize() {
    // Purpose and Method:  Check the results of sampling the Poisson and Gaussian
    //   distributions.  Calculate the mean and standard deviation.
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize" << endreq;
    

    return StatusCode::SUCCESS;
}

void TestAcdUtil::writeTileFrame(unsigned face, unsigned row, 
                                 unsigned col, bool bent) {
  idents::AcdId aid(0, face, row, col);
  idents::VolumeIdentifier vid = aid.volId(bent);

  AcdFrameUtil::AcdReferenceFrame frame =
    m_acdGeoSvc->getReferenceFrame(vid);
  (*m_log) << MSG::INFO << vid.name() << " has reference frame "
           << frame << endreq;

}
// for orient, 0 = xmeas; 1 = ymeas
void TestAcdUtil::writeRibbonFrame(unsigned face, unsigned orient, 
                                   unsigned rib, unsigned seg) {
  idents::VolumeIdentifier vid;

  vid.append(1);   // ACD
  vid.append(face);   // face
  vid.append(41);  // ribbon
  vid.append(orient);   // y-measuring
  vid.append(rib);   // ribbon number.  
  vid.append(seg);   // segment

  AcdFrameUtil::AcdReferenceFrame frame 
    =  m_acdGeoSvc->getReferenceFrame(vid);
  (*m_log) << MSG::INFO << vid.name() << " has reference frame "
           << frame << endreq;
}






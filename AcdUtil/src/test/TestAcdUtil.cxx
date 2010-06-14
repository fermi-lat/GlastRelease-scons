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

  IAcdGeometrySvc *m_acdGeoSvc;

};

static const AlgFactory<TestAcdUtil>  Factory;
const IAlgFactory& TestAcdUtilFactory = Factory;

TestAcdUtil::TestAcdUtil(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {
    
}

StatusCode TestAcdUtil::initialize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    StatusCode sc = StatusCode::SUCCESS;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    sc = service("AcdGeometrySvc", m_acdGeoSvc, true);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to find ACD Geometry service" << endreq;
        return sc;
    }

    return sc;
}


StatusCode TestAcdUtil::execute() {
    // Purpose and Method:  Exercise the methods contained in AcdDigiUtil
    //  When sampling the poisson and gaussian distributions, store the results to
    //  the appropriate vector.

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

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

        std::vector<double> dim;
        HepPoint3D center;
 
        sc = m_acdGeoSvc->getDetectorDimensions(volid, dim, center);
        if (sc.isFailure()) {
           log << MSG::WARNING << "Failed to find dimensions for " <<
               volid.name() << endreq;
           return sc;
        } 
        idents::AcdId acdId(volid);
        log << MSG::INFO << "AcdId: " << acdId.id() << " Dimensions: " 
            << dim[0] << ", " << dim[1] << ", " << dim[2] << endreq;

        HepPoint3D corner[4];
        sc = m_acdGeoSvc->getCorners(dim, center, corner);
        if (sc.isFailure()) {
           log << MSG::WARNING << "Failed to find corners for " <<
               volid.name() << endreq;
           return sc;
        } 
        unsigned int iCorner;
        for (iCorner = 0; iCorner<4; iCorner++) {
            log << "Corner " << iCorner << " " << corner[iCorner].x() 
                << ", " << corner[iCorner].y() << ", " << corner[iCorner].z()
                << endreq;
        }

    }

    log << MSG::INFO << "Number of Detectors:  " << acdDetectorCounter << endreq;

    return sc;
}


StatusCode TestAcdUtil::finalize() {
    // Purpose and Method:  Check the results of sampling the Poisson and Gaussian
    //   distributions.  Calculate the mean and standard deviation.
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize" << endreq;
    

    return StatusCode::SUCCESS;
}







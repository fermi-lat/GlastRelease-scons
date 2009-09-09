#define AcdDigi_TestAcdDigiAlg_CXX

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
#include "AcdUtil/IAcdCalibSvc.h"

#include "../../AcdDigiUtil.h"

/** @class TestAcdDigiAlg.cpp
 * @brief AcdDigiUtil test algorithm
 *
 * Exercise all of AcdDigiUtil to be sure that the methods function properly.
 *
 * $Header$
 */

class TestAcdDigiAlg : public Algorithm {

public:
  TestAcdDigiAlg(const std::string& name, ISvcLocator* pSvcLocator); 

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:

    /// local member for access to AcdDigiUtil
    AcdDigiUtil m_util;

    /// input XML file containing parameters for Digitization
    std::string	m_xmlFile;

    /// vectors to store results of sampling
    std::vector<long> m_poisson;
    std::vector<float> m_gauss;
};

static const AlgFactory<TestAcdDigiAlg>  Factory;
const IAlgFactory& TestAcdDigiAlgFactory = Factory;

TestAcdDigiAlg::TestAcdDigiAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {
    
    // Declare the properties that may be set in the job options file
    declareProperty ("xmlFile", m_xmlFile="$(ACDDIGIXMLPATH)/acdDigi.xml");
}

StatusCode TestAcdDigiAlg::initialize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    StatusCode sc = StatusCode::SUCCESS;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    IAcdGeometrySvc* acdGeomSvc(0);
    sc = service("AcdGeometrySvc", acdGeomSvc, true);
    if (sc.isSuccess() ) {
      sc = acdGeomSvc->queryInterface(IID_IAcdGeometrySvc, (void**)&acdGeomSvc);
    }
   
    if ( !sc.isSuccess() ) {
      log << MSG::WARNING << "AcdDigiAlg failed to get the AcdGeometrySvc" << endreq;
      return sc;
    } else {
      log << MSG::INFO << "Got AcdGeometrySvc" << endreq;
    }

    AcdUtil::IAcdCalibSvc* calibSvc(0);
    sc = service("AcdSimCalibSvc", calibSvc, true);
    if (sc.isSuccess() ) {
      sc = calibSvc->queryInterface(IID_IAcdCalibSvc, (void**)&calibSvc);
    }

    if ( !sc.isSuccess() ) {
      log << MSG::WARNING << "Could not get AcdSimCalibSvc"  << endreq;
      return sc;
    } else {
      log << MSG::INFO << "Got AcdSimCalibSvc" << endreq;
    } 

    m_util.initialize(*calibSvc,*acdGeomSvc,m_xmlFile);

    return sc;
}


StatusCode TestAcdDigiAlg::execute() {
    // Purpose and Method:  Exercise the methods contained in AcdDigiUtil
    //  When sampling the poisson and gaussian distributions, store the results to
    //  the appropriate vector.

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    const double energy1MeV(1.0);
    const double energy20MeV(20.0);

    idents::AcdId id(0, 1, 1, 1);

    double pe_pmt[2];
    double pmt_mips[2];
    Event::AcdDigi::Range range[2];
    unsigned short pha[2];

    sc = m_util.photoElectronsFromEnergy(id,energy1MeV,pe_pmt[0],pe_pmt[1]);
    if ( sc.isFailure() ) return sc;
    log << MSG::DEBUG << "1.0 MeV = " << pe_pmt[0] << ',' << pe_pmt[1] << " photo electrons." << endreq;
    sc = m_util.mipEquivalentLightYeild(id,pe_pmt[0],pe_pmt[1],log,pmt_mips[0],pmt_mips[0]);
    if ( sc.isFailure() ) return sc;
    log << MSG::DEBUG << "1.0 MeV = " << pmt_mips[0] << ',' << pmt_mips[0] << " mips." << endreq;
    sc = m_util.phaCounts(id,pmt_mips,false,log,range,pha);
    log << MSG::DEBUG << "1.0 MeV = " << range[0] << ':' << pha[0] << ","
	<< range[1] << ':' << pha[1] << " PHA." << endreq;
    
    sc = m_util.photoElectronsFromEnergy(id,energy20MeV,pe_pmt[0],pe_pmt[1]);
    if ( sc.isFailure() ) return sc;
    log << MSG::DEBUG << "20.0 MeV = " << pe_pmt[0] << ',' << pe_pmt[1] << " photo electrons." << endreq;
    sc = m_util.mipEquivalentLightYeild(id,pe_pmt[0],pe_pmt[1],log,pmt_mips[0],pmt_mips[0]);
    if ( sc.isFailure() ) return sc;
    log << MSG::DEBUG << "20.0 MeV = " << pmt_mips[0] << ',' << pmt_mips[0] << " mips." << endreq;
    sc = m_util.phaCounts(id,pmt_mips,false,log,range,pha);
    log << MSG::DEBUG << "20.0 MeV = " << range[0] << ':' << pha[0] << ","
	<< range[1] << ':' << pha[1] << " PHA." << endreq;

    m_poisson.push_back((long)floor(m_util.shootPoisson(5)));

    m_gauss.push_back(m_util.shootGaussian(1.0));

    return sc;
}


StatusCode TestAcdDigiAlg::finalize() {
    // Purpose and Method:  Check the results of sampling the Poisson and Gaussian
    //   distributions.  Calculate the mean and standard deviation.
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize" << endreq;
    
    float poissonMean = 0.;
    std::vector<long>::iterator poissonIt;
    for (poissonIt = m_poisson.begin(); poissonIt != m_poisson.end(); poissonIt++) {
        poissonMean += (*poissonIt);
    }
    float poissonStdDev = 0.0;
    poissonMean /= m_poisson.size();
    for (poissonIt = m_poisson.begin(); poissonIt != m_poisson.end(); poissonIt++) {
        poissonStdDev += (*poissonIt - poissonMean) * (*poissonIt - poissonMean);
    }
    poissonStdDev /= m_poisson.size();
    log << MSG::DEBUG << "Poisson with actual mean 5 returned mean: " << poissonMean
        << " stdDev: " << poissonStdDev << endreq;

    std::vector<float>::iterator gaussIt;
    float gaussMean = 0.;
    for (gaussIt = m_gauss.begin(); gaussIt != m_gauss.end(); gaussIt++) {
        gaussMean += (*gaussIt);
    }
    float gaussStdDev = 0.0;
    gaussMean /= m_gauss.size();
    for (gaussIt = m_gauss.begin(); gaussIt != m_gauss.end(); gaussIt++) {
        gaussStdDev += (*gaussIt - gaussMean) * (*gaussIt - gaussMean);
    }
    gaussStdDev /= m_gauss.size();

    log << MSG::DEBUG << "Gaussian with actual stdDev 1.0 returned mean: " << gaussMean
        << " stdDev: " << gaussStdDev << endreq;


    return StatusCode::SUCCESS;
}






